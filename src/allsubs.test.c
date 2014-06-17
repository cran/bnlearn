#include "common.h"

#define DEBUGGING() \
        if (pvalue > a) { \
          if (debuglevel > 0) { \
            Rprintf("    > node %s is independent from %s given ", \
              CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0))); \
            for (i = 0; i < cursize; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(sx, subset[i]))); \
            Rprintf("(p-value: %g).\n", pvalue); \
          }/*THEN*/ \
          UNPROTECT(3); \
          return ScalarReal(pvalue); \
        }/*THEN*/ \
        else { \
          if (debuglevel > 0) { \
            Rprintf("    > node %s is dependent on %s given ", \
              CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0))); \
            for (i = 0; i < cursize; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(sx, subset[i]))); \
            Rprintf("(p-value: %g).\n", pvalue); \
          }/*THEN*/ \
        }/*ELSE*/

#define DISCRETE_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = (int **) alloc1dpointer(nsx); \
    for (i = 0; i < nsx; i++) \
      column[i] = INTEGER(VECTOR_ELT(zz, i)); \
    /* allocate and compute the number of levels. */ \
    nlvls = alloc1dcont(nsx); \
    for (i = 0; i < nsx; i++) \
      nlvls[i] = NLEVELS(VECTOR_ELT(zz, i)); \
    /* allocate the parents' configurations. */ \
    zptr = alloc1dcont(nobs);

#define GAUSSIAN_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = (double **) alloc1dpointer(nsx + 2); \
    column[0] = REAL(xx); \
    column[1] = REAL(yy); \
    for (i = 0; i < nsx + 0; i++) \
      column[i + 2] = REAL(VECTOR_ELT(zz, i)); \
    /* allocate and compute mean values and the covariance matrix. */ \
    mean = alloc1dreal(nsx + 2); \
    c_meanvec(column, mean, nobs, nsx + 2, 0);

#define ALLOC_DISCRETE_SUBSET() \
      /* allocate and initialize the subset. */ \
      subset = Calloc(cursize, int); \
      subcol = Calloc(cursize, int *); \
      sublvls = Calloc(cursize, int); \
      /* initialize the first subset. */ \
      first_subset(subset, cursize, 0);

#define ALLOC_GAUSSIAN_SUBSET() \
      /* allocate and initialize the subset indexes array. */ \
      subset = Calloc(cursize, int); \
      /* allocate the mean values of the subset. */ \
      submean = Calloc(cursize + 2, double); \
      submean[0] = mean[0]; \
      submean[1] = mean[1]; \
      /* allocate column pointers for the subset. */ \
      subcol = Calloc(cursize + 2, double *); \
      subcol[0] = column[0]; \
      subcol[1] = column[1]; \
      /* allocate the covariance matrix and the U, D, V matrix. */ \
      cov = Calloc((cursize + 2) * (cursize + 2), double); \
      u = Calloc((cursize + 2) * (cursize + 2), double); \
      d = Calloc(cursize + 2, double); \
      vt = Calloc((cursize + 2) * (cursize + 2), double); \
      /* initialize the first subset. */ \
      first_subset(subset, cursize, 0);

#define FREE_DISCRETE_SUBSET() \
      Free(subcol); \
      Free(sublvls); \
      Free(subset);

#define FREE_GAUSSIAN_SUBSET() \
      Free(u); \
      Free(d); \
      Free(vt); \
      Free(cov); \
      Free(subset); \
      Free(submean); \
      Free(subcol);

#define PREPARE_DISCRETE_SUBSET() \
        /* prepare the variables in the current subset. */ \
        for (i = 0; i < cursize; i++) { \
          subcol[i] = column[subset[i]]; \
          sublvls[i] = nlvls[subset[i]]; \
        }/*FOR*/ \
        /* construct the parents' configurations. */ \
        cfg3(subcol, nobs, cursize, sublvls, zptr, &llz);

#define PREPARE_GAUSSIAN_SUBSET() \
        /* prepare the variables in the current subset. */ \
        for (i = 0; i < cursize; i++) { \
          subcol[i + 2] = column[subset[i] + 2]; \
          submean[i + 2] = mean[subset[i] + 2]; \
        }/*FOR*/

SEXP allsubs_test(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP min, SEXP max, SEXP debug) {

int minsize = INT(min), maxsize = INT(max), debuglevel = isTRUE(debug);
int i = 0, *subset = NULL, cursize = 0, nsx = length(sx), nobs = 0;
double statistic = 0, pvalue = 0, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
void *memp = NULL;
SEXP xx, yy, zz;

  /* call utest to deal with zero-length conditioning subsets. */
  if (minsize == 0) {

    pvalue = NUM(utest(x, y, data, test, B, alpha, TRUESEXP));

    if (pvalue > a) {

      if (debuglevel > 0)
        Rprintf("    > node %s is marginally independent from %s (p-value: %g).\n",
          CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)), pvalue);

      return ScalarReal(pvalue);

    }/*THEN*/
    else {

      if (debuglevel > 0)
        Rprintf("    > node %s is marginally dependent on %s (p-value: %g).\n",
          CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)), pvalue);

    }/*ELSE*/

  }/*THEN*/

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, TRUE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, FALSE));
  nobs = length(xx);

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL;
    int llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
    int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
    int adj = (strcmp(t, "mi-adf") == 0);
    double df = 0;

    DISCRETE_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_DISCRETE_SUBSET();

      /* iterate over subsets. */
      do {

        PREPARE_DISCRETE_SUBSET();

        memp = vmaxget();

        statistic = c_cmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df, adj);
        pvalue = pchisq(2 * nobs * statistic, df, FALSE, FALSE);

        vmaxset(memp);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_DISCRETE_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL;
    int llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
    int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
    double df = 0;

    DISCRETE_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_DISCRETE_SUBSET();

      /* iterate over subsets. */
      do {

        PREPARE_DISCRETE_SUBSET();

        memp = vmaxget();

        statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df);
        pvalue = pchisq(2 * nobs * statistic, df, FALSE, FALSE);

        vmaxset(memp);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_DISCRETE_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL;
    int llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
    int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
    int adj = (strcmp(t, "x2-adf") == 0);
    double df = 0;

    DISCRETE_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_DISCRETE_SUBSET();

      /* iterate over subsets. */
      do {

        PREPARE_DISCRETE_SUBSET();

        memp = vmaxget();

        statistic = c_cx2(xptr, llx, yptr, lly, zptr, llz, nobs, &df, adj);
        pvalue = pchisq(statistic, df, FALSE, FALSE);

        vmaxset(memp);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_DISCRETE_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL;
    int llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
    int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;

    DISCRETE_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_DISCRETE_SUBSET();

      /* iterate over subsets. */
      do {

        PREPARE_DISCRETE_SUBSET();

        memp = vmaxget();

        statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
        pvalue = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

        vmaxset(memp);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_DISCRETE_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* Pearson's linear correlation test. */
    double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
    double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;

    GAUSSIAN_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      /* check that we have enough degrees of freedom. */
      if (nobs - cursize - 2 < 1)
        error("trying to do a conditional independence test with zero degrees of freedom.");

      ALLOC_GAUSSIAN_SUBSET();

      do {

        PREPARE_GAUSSIAN_SUBSET();

        /* compute the covariance matrix. */
        c_covmat(subcol, submean, cursize + 2, nobs, cov, 0);
        statistic = c_fast_pcor(cov, u, d, vt, cursize + 2, TRUE);
        statistic = fabs(statistic * sqrt(nobs - cursize - 2) / 
                            sqrt(1 - statistic * statistic));
        pvalue = 2 * pt(statistic, nobs - cursize - 2, FALSE, FALSE);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_GAUSSIAN_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
    double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;

    GAUSSIAN_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      /* check that we have enough degrees of freedom. */
      if (nobs - cursize - 2 < 1)
        error("trying to do a conditional independence test with zero degrees of freedom.");

      ALLOC_GAUSSIAN_SUBSET();

      do {

        PREPARE_GAUSSIAN_SUBSET();

        /* compute the covariance matrix. */
        c_covmat(subcol, submean, cursize + 2, nobs, cov, 0);
        statistic = c_fast_pcor(cov, u, d, vt, cursize + 2, TRUE);
        statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                      sqrt((double)nobs - cursize - 3);
        pvalue = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_GAUSSIAN_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
    double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;

    GAUSSIAN_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_GAUSSIAN_SUBSET();

      do {

        PREPARE_GAUSSIAN_SUBSET();

        /* compute the covariance matrix. */
        c_covmat(subcol, submean, cursize + 2, nobs, cov, 0);
        statistic = c_fast_pcor(cov, u, d, vt, cursize + 2, TRUE);
        statistic = - nobs * log(1 - statistic * statistic);
        pvalue = pchisq(statistic, 1, FALSE, FALSE); 

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_GAUSSIAN_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
    double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;

    GAUSSIAN_CACHE();

    for (cursize = 1; cursize <= maxsize; cursize++) {

      ALLOC_GAUSSIAN_SUBSET();

      do {

        PREPARE_GAUSSIAN_SUBSET();

        /* compute the covariance matrix. */
        c_cov_lambda(subcol, submean, cursize + 2, nobs, cov);
        statistic = c_fast_pcor(cov, u, d, vt, cursize + 2, TRUE);
        statistic = - nobs * log(1 - statistic * statistic);
        pvalue = pchisq(statistic, 1, FALSE, FALSE);

        DEBUGGING();

      } while (next_subset(subset, cursize, nsx, 0) != NULL);

      FREE_GAUSSIAN_SUBSET();

    }/*FOR*/

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0;
    double a2 = (strncmp(t, "smc-", 4) == 0) ? a : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type)) {

      int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL;
      int llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
      int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
      double df = 0;
      
      DISCRETE_CACHE();

      for (cursize = 1; cursize <= maxsize; cursize++) {

        ALLOC_DISCRETE_SUBSET();

        /* iterate over subsets. */
        do {

          PREPARE_DISCRETE_SUBSET();

          memp = vmaxget();

          c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, INT(B), &statistic,
            &pvalue, a2, type, &df);

          vmaxset(memp);

          DEBUGGING();

        } while (next_subset(subset, cursize, nsx, 0) != NULL);

        FREE_DISCRETE_SUBSET();

      }/*FOR*/

    }/*THEN*/
    else {

      /* Gaussian test. */
      double **column = NULL, **subcol= NULL;

      /* allocate and initialize an array of pointers for the variables. */
      column = (double **) alloc1dpointer(nsx + 2);
      column[0] = REAL(xx);
      column[1] = REAL(yy);
      for (i = 0; i < nsx + 0; i++)
        column[i + 2] = REAL(VECTOR_ELT(zz, i));

      for (cursize = 1; cursize <= maxsize; cursize++) {

        /* allocate and initialize the subset indexes array. */
        subset = Calloc(cursize, int);
        /* */
        subcol = Calloc(cursize + 2, double *);
        subcol[0] = column[0];
        subcol[1] = column[1];
        /* initialize the first subset. */
        first_subset(subset, cursize, 1);

        /* iterate over subsets. */
        do {

          /* prepare the variables in the current subset. */
          for (i = 0; i < cursize; i++)
            subcol[i + 2] = column[subset[i] + 1];

          memp = vmaxget();

          c_gauss_cmcarlo(subcol, cursize + 2, nobs, INT(B), &statistic, &pvalue,
            a2, type);

          vmaxset(memp);

          DEBUGGING();

        } while (next_subset(subset, cursize, nsx, 1) != NULL);

        Free(subset);
        Free(subcol);

      }/*FOR*/

    }/*ELSE*/

  }/*THEN*/

  UNPROTECT(3);

  return ScalarReal(pvalue);

}/*ALLSUBS_TEST*/

