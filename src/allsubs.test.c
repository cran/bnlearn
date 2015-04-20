#include "include/rcore.h"
#include "include/allocations.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/dataframe.h"
#include "include/globals.h"
#include "include/covariance.h"

#define DEBUGGING(UNWIND) \
        if (pvalue > a) { \
          if (debuglevel > 0) { \
            Rprintf("    > node %s is independent from %s %s", \
              CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)), \
              (nf > 0 || cursize > 0) ? "given " : ""); \
            for (i = 0; i < nf; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(fixed, i))); \
            for (i = nf; i < cursize + nf; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(sx, subset[i] - nf))); \
            Rprintf("(p-value: %g).\n", pvalue); \
          }/*THEN*/ \
          UNWIND; \
          return mkRealVec(3, pvalue, min_pvalue, max_pvalue); \
        }/*THEN*/ \
        else { \
          if (debuglevel > 0) { \
            Rprintf("    > node %s is dependent on %s %s", \
              CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)), \
              (nf > 0 || cursize > 0) ? "given " : ""); \
            for (i = 0; i < nf; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(fixed, i))); \
            for (i = nf; i < cursize + nf; i++) \
              Rprintf("%s ", CHAR(STRING_ELT(sx, subset[i] - nf))); \
            Rprintf("(p-value: %g).\n", pvalue); \
          }/*THEN*/ \
        }/*ELSE*/

#define DISCRETE_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = (int **) alloc1dpointer(nsx + nf); \
    for (i = 0; i < nf; i++) \
      column[i] = INTEGER(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx; i++) \
      column[i + nf] = INTEGER(VECTOR_ELT(zz, i)); \
    /* allocate and compute the number of levels. */ \
    nlvls = alloc1dcont(nsx + nf); \
    for (i = 0; i < nf; i++) \
      nlvls[i] = NLEVELS(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx; i++) \
      nlvls[i + nf] = NLEVELS(VECTOR_ELT(zz, i)); \
    /* allocate the parents' configurations. */ \
    zptr = alloc1dcont(nobs);

#define GAUSSIAN_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = (double **) alloc1dpointer(nsx + nf + 2); \
    column[0] = REAL(xx); \
    column[1] = REAL(yy); \
    for (i = 0; i < nf; i++) \
      column[i + 2] = REAL(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx + 0; i++) \
      column[i + nf + 2] = REAL(VECTOR_ELT(zz, i)); \
    /* allocate and compute mean values and the covariance matrix. */ \
    mean = alloc1dreal(nsx + nf + 2); \
    c_meanvec(column, mean, nobs, nsx + nf + 2, 0);

#define ALLOC_DISCRETE_SUBSET() \
      /* allocate and initialize the subset. */ \
      subset = Calloc(cursize + nf, int); \
      subcol = Calloc(cursize + nf, int *); \
      sublvls = Calloc(cursize + nf, int); \
      /* initialize the first subset. */ \
      first_subset(subset + nf, cursize, nf); \
      for (i = 0; i < nf; i++) \
        subset[i] = i;

#define ALLOC_GAUSSIAN_SUBSET() \
      /* allocate and initialize the subset indexes array. */ \
      subset = Calloc(cursize + nf, int); \
      /* allocate the mean values of the subset. */ \
      submean = Calloc(cursize + nf + 2, double); \
      submean[0] = mean[0]; \
      submean[1] = mean[1]; \
      for (i = 0; i < nf; i++) \
        submean[i + 2] = mean[i + 2]; \
      /* allocate column pointers for the subset. */ \
      subcol = Calloc(cursize + nf + 2, double *); \
      subcol[0] = column[0]; \
      subcol[1] = column[1]; \
      for (i = 0; i < nf; i++) \
        subcol[i + 2] = column[i + 2]; \
      /* allocate the covariance matrix and the U, D, V matrix. */ \
      cov = Calloc((cursize + nf + 2) * (cursize + nf + 2), double); \
      u = Calloc((cursize + nf + 2) * (cursize + nf + 2), double); \
      d = Calloc(cursize + nf + 2, double); \
      vt = Calloc((cursize + nf + 2) * (cursize + nf + 2), double); \
      /* initialize the first subset. */ \
      first_subset(subset + nf, cursize, nf); \
      for (i = 0; i < nf; i++) \
        subset[i] = i;

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
        for (i = 0; i < cursize + nf; i++) { \
          subcol[i] = column[subset[i]]; \
          sublvls[i] = nlvls[subset[i]]; \
        }/*FOR*/ \
        /* construct the parents' configurations. */ \
        c_fast_config(subcol, nobs, cursize + nf, sublvls, zptr, &llz, 1);

#define PREPARE_GAUSSIAN_SUBSET() \
        /* prepare the variables in the current subset. */ \
        for (i = 0; i < cursize + nf; i++) { \
          subcol[i + 2] = column[subset[i] + 2]; \
          submean[i + 2] = mean[subset[i] + 2]; \
        }/*FOR*/

#define PVALUE(testfun) \
        pvalue = testfun; \
        min_pvalue = pvalue < min_pvalue ? pvalue : min_pvalue; \
        max_pvalue = pvalue > max_pvalue ? pvalue : max_pvalue;

/* mutual information test, with and without df adjustments. */
static inline SEXP ast_mi(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int adj, int a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
void *memp = NULL;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      memp = vmaxget();
      statistic = c_cmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df, adj);
      PVALUE(pchisq(2 * nobs * statistic, df, FALSE, FALSE));
      vmaxset(memp);
      DEBUGGING(FREE_DISCRETE_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MI*/

/* shrinkage mutual information test. */
static inline SEXP ast_mish(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;void *memp = NULL;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      memp = vmaxget();
      statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df);
      PVALUE(pchisq(2 * nobs * statistic, df, FALSE, FALSE));
      vmaxset(memp);
      DEBUGGING(FREE_DISCRETE_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MISH*/

/* Pearson's X^2 test, with and without df adjustments. */
static inline SEXP ast_x2(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int adj, int a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
void *memp = NULL;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      memp = vmaxget();
      statistic = c_cx2(xptr, llx, yptr, lly, zptr, llz, nobs, &df, adj);
      PVALUE(pchisq(statistic, df, FALSE, FALSE));
      vmaxset(memp);
      DEBUGGING(FREE_DISCRETE_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_X2*/

/* Jonckheere-Terpstra test. */
static inline SEXP ast_jt(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
void *memp = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      memp = vmaxget();
      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
      PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));
      vmaxset(memp);
      DEBUGGING(FREE_DISCRETE_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_JT*/

/* Pearson's linear correlation test. */
static inline SEXP ast_cor(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* check that we have enough degrees of freedom. */
    if (nobs - cursize - nf - 2 < 1)
      error("trying to do a conditional independence test with zero degrees of freedom.");

    ALLOC_GAUSSIAN_SUBSET();

    do {

      PREPARE_GAUSSIAN_SUBSET();
      /* compute the covariance matrix. */
      c_covmat(subcol, submean, cursize + nf + 2, nobs, cov, 0);
      statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
      statistic = fabs(statistic * sqrt(nobs - cursize - nf - 2) /
                          sqrt(1 - statistic * statistic));
      PVALUE(2 * pt(statistic, nobs - cursize - nf - 2, FALSE, FALSE));
      DEBUGGING(FREE_GAUSSIAN_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_GAUSSIAN_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_COR*/

/* Fisher's Z test. */
static inline SEXP ast_zf(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* check that we have enough degrees of freedom. */
    if (nobs - cursize - nf - 2 < 1)
      error("trying to do a conditional independence test with zero degrees of freedom.");

    ALLOC_GAUSSIAN_SUBSET();

    do {

      PREPARE_GAUSSIAN_SUBSET();
      /* compute the covariance matrix. */
      c_covmat(subcol, submean, cursize + nf + 2, nobs, cov, 0);
      statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
      statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                    sqrt((double)nobs - cursize - nf - 3);
      PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));
      DEBUGGING(FREE_GAUSSIAN_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_GAUSSIAN_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_ZF*/

/* Gaussian mutual information test. */
static inline SEXP ast_mig(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_GAUSSIAN_SUBSET();

    do {

      PREPARE_GAUSSIAN_SUBSET();
      /* compute the covariance matrix. */
      c_covmat(subcol, submean, cursize + nf + 2, nobs, cov, 0);
      statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
      statistic = - nobs * log(1 - statistic * statistic);
      PVALUE(pchisq(statistic, 1, FALSE, FALSE));
      DEBUGGING(FREE_GAUSSIAN_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_GAUSSIAN_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MIG*/

/* shrinkage Gaussian mutual information test. */
static inline SEXP ast_migsh(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_GAUSSIAN_SUBSET();

    do {

      PREPARE_GAUSSIAN_SUBSET();
      /* compute the covariance matrix. */
      c_cov_lambda(subcol, submean, cursize + nf + 2, nobs, cov);
      statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
      statistic = - nobs * log(1 - statistic * statistic);
      PVALUE(pchisq(statistic, 1, FALSE, FALSE));
      DEBUGGING(FREE_GAUSSIAN_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_GAUSSIAN_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MIGSH*/

/* conditional linear Gaussian test. */
static inline SEXP ast_micg(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int debuglevel) {

int cursize = 0, ndp = 0, ngp = 0, xtype = TYPEOF(xx), ytype = TYPEOF(yy);
int i = 0, j = 0, k = 0, *subset = NULL, *nlvls = NULL, **dp = NULL;
int *zptr = NULL, llx = 0, lly = 0, llz = 0, *dlvls = NULL;
double **gp = NULL, statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
void *xptr = 0, *yptr = 0, **columns = NULL;
SEXP temp;

  /* scan the parents for type and number of levels. */
  columns = Calloc(nsx + nf, void *);
  nlvls = Calloc(nsx + nf, int);

  for (i = 0; i < nf; i++) {

    temp = VECTOR_ELT(ff, i);

    if (TYPEOF(temp) == INTSXP) {

      columns[i] = INTEGER(temp);
      nlvls[i] = NLEVELS(temp);
      ndp++;

    }/*THEN*/
    else {

      columns[i] = REAL(temp);
      ngp++;

    }/*ELSE*/

  }/*FOR*/

  for (i = 0; i < nsx; i++) {

    temp = VECTOR_ELT(zz, i);

    if (TYPEOF(temp) == INTSXP) {

      columns[i + nf] = INTEGER(temp);
      nlvls[i + nf] = NLEVELS(temp);
      ndp++;

    }/*THEN*/
    else {

      columns[i + nf] = REAL(temp);
      ngp++;

    }/*ELSE*/

  }/*FOR*/

  /* if both variables are continuous and all conditioning variables are
   * continuous,  the test reverts back to a Gaussian mutual information
   * test. */
  if ((xtype == INTSXP) && (ytype == INTSXP)) {

    xptr = INTEGER(xx);
    llx = NLEVELS(xx);
    yptr = INTEGER(yy);
    lly = NLEVELS(yy);

  }/*THEN*/
  if ((xtype == REALSXP) && (ytype == REALSXP)) {

    xptr = REAL(xx);
    yptr = REAL(yy);

  }/*THEN*/
  else if ((xtype == REALSXP) && (ytype == INTSXP)) {

    xptr = REAL(xx);
    yptr = INTEGER(yy);
    lly = NLEVELS(yy);

  }/*THEN*/
  else if ((xtype == INTSXP) && (ytype == REALSXP)) {

    yptr = INTEGER(xx);
    lly = NLEVELS(xx);
    xptr = REAL(yy);

  }/*THEN*/

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset. */
    subset = Calloc(cursize + nf, int);
    /* initialize the first subset. */
    first_subset(subset + nf, cursize, nf);
    for (i = 0; i < nf; i++)
      subset[i] = i;

    do {

      /* count how many discrete and continuous parents are in the subset. */
      for (i = 0, ndp = 0; i < cursize + nf; i++)
        ndp += (nlvls[subset[i]] > 0);
      ngp = cursize + nf - ndp;

      /* split parents by type. */
      gp = Calloc(ngp + 1, double *);
      dp = Calloc(ndp + 1, int *);
      dlvls = Calloc(ndp + 1, int);
      for (i = 0, j = 0, k = 0; i < cursize + nf; i++)
        if (nlvls[subset[i]] > 0) {

          dp[1 + j] = columns[subset[i]];
          dlvls[1 + j++] = nlvls[subset[i]];

        }/*THEN*/
        else {

          gp[1 + k++] = columns[subset[i]];

        }/*ELSE*/

      /* if there are discrete conditioning variables, compute their
       * configurations. */
      if (ndp > 0) {

        zptr = Calloc(nobs, int);
        c_fast_config(dp + 1, nobs, ndp, dlvls + 1, zptr, &llz, 1);

      }/*THEN*/

      if ((ytype == INTSXP) && (xtype == INTSXP)) {

        /* check whether the conditioning set is valid. */
        if (ngp > 0) {

          /* need to reverse conditioning to actually compute the test. */
          statistic = 2 * nobs * c_cmicg_unroll(xptr, llx, yptr, lly, zptr, llz,
                                   gp + 1, ngp, &df, nobs);

        }/*THEN*/
        else {

          /* if both nodes are discrete, the test reverts back to a discrete
           * mutual information test. */
          statistic = 2 * nobs * c_cmi(xptr, llx, yptr, lly, zptr, llz,
                                   nobs, &df, FALSE);

        }/*ELSE*/

      }/*THEN*/
      else if ((ytype == REALSXP) && (xtype == REALSXP)) {

        gp[0] = xptr;
        statistic = 2 * nobs * c_cmicg(yptr, gp, ngp + 1, NULL, 0, zptr, llz,
                               dlvls, nobs);

        /* one regression coefficient for each conditioning level is added;
         * if all conditioning variables are continuous that's just one global
         * regression coefficient. */
        df = (llz == 0) ? 1 : llz;

      }/*THEN*/
      else { /* one variable is discrete, the other is continuous. */

        dp[0] = yptr;
        dlvls[0] = lly;
        statistic = 2 * nobs * c_cmicg(xptr, gp + 1, ngp, dp, ndp + 1, zptr,
                                 llz, dlvls, nobs);

        /* for each additional configuration of the discrete conditioning
         * variables plus the discrete yptr, one whole set of regression
         * coefficients (plus the intercept) is added. */
        df = (lly - 1) * ((llz == 0) ? 1 : llz)  * (ngp + 1);

      }/*THEN*/

      Free(gp);
      Free(dp);
      Free(dlvls);
      Free(zptr);

      PVALUE(pchisq(statistic, df, FALSE, FALSE));
      DEBUGGING(Free(subset); Free(columns); Free(nlvls); Free(subset));

    } while (next_subset(subset + nf, cursize, nsx, nf));

    Free(subset);

  }/*FOR*/

  Free(columns);
  Free(nlvls);

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MICG*/

/* discrete permutation tests. */
static inline SEXP ast_dperm(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int type, int nperms, int threshold, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
void *memp = NULL;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      memp = vmaxget();
      c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, nperms, &statistic,
        &pvalue, threshold, type, &df);
      PVALUE(pvalue);
      vmaxset(memp);
      DEBUGGING(FREE_DISCRETE_SUBSET());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_DPERM*/

/* continuous permutation tests. */
static inline SEXP ast_gperm(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    int a, int type, int nperms, int threshold, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, **subcol= NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
void *memp = NULL;

  /* allocate and initialize an array of pointers for the variables. */
  column = (double **) alloc1dpointer(nsx + nf + 2);
  column[0] = REAL(xx);
  column[1] = REAL(yy);
  for (i = 0; i < nf; i++)
    column[i + 2] = REAL(VECTOR_ELT(ff, i));
  for (i = 0; i < nsx + 0; i++)
    column[i + nf + 2] = REAL(VECTOR_ELT(zz, i));

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset indexes array. */
    subset = Calloc(cursize + nf, int);
    /* allocate and initialize the subset. */
    subcol = Calloc(cursize + nf + 2, double *);
    subcol[0] = column[0];
    subcol[1] = column[1];
    for (i = 0; i < nf; i++)
      subcol[i + 2] = REAL(VECTOR_ELT(ff, i));
    /* initialize the first subset. */
    first_subset(subset + nf, cursize, nf);
    for (i = 0; i < nf; i++)
      subset[i] = i;

    /* iterate over subsets. */
    do {

      /* prepare the variables in the current subset. */
      for (i = 0; i < cursize + nf; i++)
        subcol[i + 2] = column[subset[i] + 2];
      memp = vmaxget();
      c_gauss_cmcarlo(subcol, cursize + nf + 2, nobs, nperms, &statistic, &pvalue,
        threshold, type);
      PVALUE(pvalue);
      vmaxset(memp);
      DEBUGGING(Free(subset); Free(subcol));

    } while (next_subset(subset + nf, cursize, nsx, nf));

    Free(subset);
    Free(subcol);

  }/*FOR*/

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_GPERM*/

SEXP allsubs_test(SEXP x, SEXP y, SEXP sx, SEXP fixed, SEXP data, SEXP test,
    SEXP B, SEXP alpha, SEXP min, SEXP max, SEXP debug) {

int minsize = INT(min), maxsize = INT(max), debuglevel = isTRUE(debug), nobs = 0;
int i = 0, *subset = NULL, cursize = 0, nsx = length(sx), nf = length(fixed);
double pvalue = 0, min_pvalue = 1, max_pvalue = 0, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
SEXP xx, yy, zz, ff = R_NilValue, res;

  /* call indep_test to deal with zero-length conditioning subsets. */
  if (minsize == 0) {

    PVALUE(NUM(indep_test(x, y, fixed, data, test, B, alpha, TRUESEXP)));
    DEBUGGING();

  }/*THEN*/

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, TRUE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, FALSE));
  if (nf > 0)
    PROTECT(ff = c_dataframe_column(data, fixed, FALSE, FALSE));
  nobs = length(xx);

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    res = ast_mi(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, (strcmp(t, "mi-adf") == 0), NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    res = ast_mish(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    res = ast_x2(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, (strcmp(t, "x2-adf") == 0), NUM(alpha), debuglevel);


  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    res = ast_jt(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* Pearson's linear correlation test. */
    res = ast_cor(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    res = ast_zf(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    res = ast_mig(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    res = ast_migsh(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if (strcmp(t, "mi-cg") == 0) {

    /* conditional linear Gaussian test. */
    res = ast_micg(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, NUM(alpha), debuglevel);

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0, nperms = INT(B);
    double threshold = (strncmp(t, "smc-", 4) == 0) ? a : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type))
      res = ast_dperm(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
              maxsize, NUM(alpha), type, nperms, threshold, debuglevel);
    else
      res = ast_gperm(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
              maxsize, NUM(alpha), type, nperms, threshold, debuglevel);

  }/*THEN*/
  else {

    error("unknown test statistic '%s'.", t);

  }/*ELSE*/

  UNPROTECT(3 + (nf > 0));

  return res;

}/*ALLSUBS_TEST*/

