#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/dataframe.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"

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
    column = (int **) Calloc1D(nsx + nf, sizeof(int *)); \
    for (i = 0; i < nf; i++) \
      column[i] = INTEGER(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx; i++) \
      column[i + nf] = INTEGER(VECTOR_ELT(zz, i)); \
    /* allocate and compute the number of levels. */ \
    nlvls = Calloc1D(nsx + nf, sizeof(int)); \
    for (i = 0; i < nf; i++) \
      nlvls[i] = NLEVELS(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx; i++) \
      nlvls[i + nf] = NLEVELS(VECTOR_ELT(zz, i)); \
    /* allocate the parents' configurations. */ \
    zptr = Calloc1D(nobs, sizeof(int));

#define FREE_DISCRETE_CACHE() \
  Free1D(column); \
  Free1D(nlvls); \
  Free1D(zptr);

#define GAUSSIAN_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = (double **) Calloc1D(nsx + nf + 2, sizeof(double *)); \
    column[0] = REAL(xx); \
    column[1] = REAL(yy); \
    for (i = 0; i < nf; i++) \
      column[i + 2] = REAL(VECTOR_ELT(ff, i)); \
    for (i = 0; i < nsx + 0; i++) \
      column[i + nf + 2] = REAL(VECTOR_ELT(zz, i));

#define FREE_GAUSSIAN_CACHE() \
  Free1D(column);

#define GAUSSIAN_ALLOC() \
    /* allocate and compute mean values and the covariance matrix. */ \
    mean = Calloc1D(nsx + nf + 2, sizeof(double)); \
    c_meanvec(column, mean, nobs, nsx + nf + 2, 0);

#define FREE_GAUSSIAN_ALLOC() \
  Free1D(mean);

#define ALLOC_DISCRETE_SUBSET() \
      /* allocate and initialize the subset. */ \
      subset = Calloc1D(cursize + nf, sizeof(int)); \
      subcol = Calloc1D(cursize + nf, sizeof(int *)); \
      sublvls = Calloc1D(cursize + nf, sizeof(int)); \
      /* initialize the first subset. */ \
      first_subset(subset + nf, cursize, nf); \
      for (i = 0; i < nf; i++) \
        subset[i] = i;

#define ALLOC_GAUSSIAN_SUBSET() \
      /* allocate and initialize the subset indexes array. */ \
      subset = Calloc1D(cursize + nf, sizeof(int)); \
      /* allocate the mean values of the subset. */ \
      submean = Calloc1D(cursize + nf + 2, sizeof(double)); \
      submean[0] = mean[0]; \
      submean[1] = mean[1]; \
      for (i = 0; i < nf; i++) \
        submean[i + 2] = mean[i + 2]; \
      /* allocate column pointers for the subset. */ \
      subcol = Calloc1D(cursize + nf + 2, sizeof(double *)); \
      subcol[0] = column[0]; \
      subcol[1] = column[1]; \
      for (i = 0; i < nf; i++) \
        subcol[i + 2] = column[i + 2]; \
      /* allocate the covariance matrix and the U, D, V matrix. */ \
      cov = Calloc1D((cursize + nf + 2) * (cursize + nf + 2), sizeof(double)); \
      c_udvt(&u, &d, &vt, (cursize + nf + 2)); \
      /* initialize the first subset. */ \
      first_subset(subset + nf, cursize, nf); \
      for (i = 0; i < nf; i++) \
        subset[i] = i;

#define FREE_DISCRETE_SUBSET() \
      Free1D(subcol); \
      Free1D(sublvls); \
      Free1D(subset);

#define FREE_GAUSSIAN_SUBSET() \
      Free1D(u); \
      Free1D(d); \
      Free1D(vt); \
      Free1D(cov); \
      Free1D(subset); \
      Free1D(submean); \
      Free1D(subcol);

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

/* parametric tests for discrete variables. */
static SEXP ast_discrete(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    test_e test, double a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();

      if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

        /* mutual information and Pearson's X^2 asymptotic tests. */
        statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, nobs, &df, test);
        if ((test == MI) || (test == MI_ADF))
          statistic = 2 * nobs * statistic;
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_SH) {

        /* shrinkage mutual information test. */
        statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df);
        PVALUE(pchisq(2 * nobs * statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == JT) {

        /* Jonckheere-Terpstra test. */
        statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
        PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));

      }/*THEN*/

      DEBUGGING(FREE_DISCRETE_SUBSET(); FREE_DISCRETE_CACHE());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  FREE_DISCRETE_CACHE();

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_DISCRETE*/

/* parametric tests for Gaussian variables. */
static SEXP ast_gaustests(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    double a, int debuglevel, test_e test) {

int i = 0, cursize = 0, *subset = NULL, df = 0;
double **column = NULL, *mean = NULL, **subcol = NULL, *submean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL, lambda = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();
  GAUSSIAN_ALLOC();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* compute the degrees of freedom for correlation and mutual information. */
    if (test == COR)
      df = nobs - (cursize + nf + 2);
    else if ((test == MI_G) || (test == MI_G_SH))
      df = 1;

    if (((test == COR) && (df < 1)) ||
        ((test == ZF) && (nobs - (cursize + nf + 2) < 2))) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      PVALUE(1);
      DEBUGGING(FREE_GAUSSIAN_CACHE(); FREE_GAUSSIAN_ALLOC());

    }/*THEN*/

    ALLOC_GAUSSIAN_SUBSET();

    do {

      PREPARE_GAUSSIAN_SUBSET();
      /* compute the covariance matrix. */
      c_covmat(subcol, submean, cursize + nf + 2, nobs, cov, 0);

      if (test == COR) {

        statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
        statistic = cor_t_trans(statistic, (double)df);
        PVALUE(2 * pt(fabs(statistic), df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_G) {

        statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
        statistic = 2 * nobs * cor_mi_trans(statistic);
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_G_SH) {

        lambda = covmat_lambda(subcol, submean, cov, nobs, cursize + nf + 2);
        covmat_shrink(cov, cursize + nf + 2, lambda);
        statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
        statistic = 2 * nobs * cor_mi_trans(statistic);
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == ZF) {

        statistic = c_fast_pcor(cov, u, d, vt, cursize + nf + 2, TRUE);
        statistic = cor_zf_trans(statistic, (double)nobs - (cursize + nf + 2));
        PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));

      }/*THEN*/

      DEBUGGING(FREE_GAUSSIAN_SUBSET(); FREE_GAUSSIAN_CACHE(); FREE_GAUSSIAN_ALLOC());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_GAUSSIAN_SUBSET();

  }/*FOR*/

  FREE_GAUSSIAN_CACHE();
  FREE_GAUSSIAN_ALLOC();

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_GAUSTESTS*/

/* conditional linear Gaussian test. */
static SEXP ast_micg(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    double a, int debuglevel) {

int cursize = 0, ndp = 0, ngp = 0, xtype = TYPEOF(xx), ytype = TYPEOF(yy);
int i = 0, j = 0, k = 0, *subset = NULL, *nlvls = NULL, **dp = NULL;
int *zptr = NULL, llx = 0, lly = 0, llz = 0, *dlvls = NULL;
double **gp = NULL, statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
void *xptr = 0, *yptr = 0, **columns = NULL;

  /* scan the parents for type and number of levels. */
  columns = Calloc1D(nsx + nf, sizeof(void *));
  nlvls = Calloc1D(nsx + nf, sizeof(int));
  df2micg(ff, columns, nlvls, &ndp, &ngp);
  df2micg(zz, columns + nf, nlvls + nf, &ndp, &ngp);

  /* if both variables are continuous and all conditioning variables are
   * continuous, the test reverts back to a Gaussian mutual information test. */
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
    subset = Calloc1D(cursize + nf, sizeof(int));
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
      gp = Calloc1D(ngp + 1, sizeof(double *));
      dp = Calloc1D(ndp + 1, sizeof(int *));
      dlvls = Calloc1D(ndp + 1, sizeof(int));
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

        zptr = Calloc1D(nobs, sizeof(int));
        c_fast_config(dp + 1, nobs, ndp, dlvls + 1, zptr, &llz, 1);

      }/*THEN*/

      if ((ytype == INTSXP) && (xtype == INTSXP)) {

        /* check whether the conditioning set is valid. */
        if (ngp > 0) {

          /* need to reverse conditioning to actually compute the test. */
          statistic = 2 * nobs * nobs *
                        c_cmicg_unroll(xptr, llx, yptr, lly, zptr, llz,
                                   gp + 1, ngp, &df, nobs);

        }/*THEN*/
        else {

          /* if both nodes are discrete, the test reverts back to a discrete
           * mutual information test. */
          statistic = 2 * nobs * c_cchisqtest(xptr, llx, yptr, lly, zptr, llz,
                                   nobs, &df, MI);

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

      }/*ELSE*/

      Free1D(gp);
      Free1D(dp);
      Free1D(dlvls);
      Free1D(zptr);

      PVALUE(pchisq(statistic, df, FALSE, FALSE));
      DEBUGGING(Free1D(subset); Free1D(columns); Free1D(nlvls); Free1D(subset));

    } while (next_subset(subset + nf, cursize, nsx, nf));

    Free1D(subset);

  }/*FOR*/

  Free1D(columns);
  Free1D(nlvls);

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_MICG*/

/* discrete permutation tests. */
static SEXP ast_dperm(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    double a, test_e type, int nperms, double threshold, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
int **column = NULL, *nlvls = NULL, **subcol = NULL, *sublvls = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;

  DISCRETE_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    ALLOC_DISCRETE_SUBSET();

    /* iterate over subsets. */
    do {

      PREPARE_DISCRETE_SUBSET();
      c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, nperms, &statistic,
        &pvalue, threshold, type, &df);
      PVALUE(pvalue);
      DEBUGGING(FREE_DISCRETE_SUBSET(); FREE_DISCRETE_CACHE());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    FREE_DISCRETE_SUBSET();

  }/*FOR*/

  FREE_DISCRETE_CACHE();

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_DPERM*/

/* continuous permutation tests. */
static SEXP ast_gperm(SEXP xx, SEXP yy, SEXP zz, SEXP ff, SEXP x, SEXP y,
    SEXP sx, SEXP fixed, int nobs, int nsx, int nf, int minsize, int maxsize,
    double a, test_e type, int nperms, double threshold, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double **column = NULL, **subcol= NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;

  GAUSSIAN_CACHE();

  for (cursize = 1; cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset indexes array. */
    subset = Calloc1D(cursize + nf, sizeof(int));
    /* allocate and initialize the subset. */
    subcol =  Calloc1D(cursize + nf + 2, sizeof(double *));
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
      c_gauss_cmcarlo(subcol, cursize + nf + 2, nobs, nperms, &statistic, &pvalue,
        threshold, type);
      PVALUE(pvalue);
      DEBUGGING(Free1D(subset); Free1D(subcol); FREE_GAUSSIAN_CACHE());

    } while (next_subset(subset + nf, cursize, nsx, nf));

    Free1D(subset);
    Free1D(subcol);

  }/*FOR*/

  FREE_GAUSSIAN_CACHE();

  return mkRealVec(3, pvalue, min_pvalue, max_pvalue);

}/*AST_GPERM*/

SEXP allsubs_test(SEXP x, SEXP y, SEXP sx, SEXP fixed, SEXP data, SEXP test,
    SEXP B, SEXP alpha, SEXP min, SEXP max, SEXP debug) {

int minsize = INT(min), maxsize = INT(max), debuglevel = isTRUE(debug), nobs = 0;
int i = 0, *subset = NULL, cursize = 0, nsx = length(sx), nf = length(fixed);
double pvalue = 0, min_pvalue = 1, max_pvalue = 0, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
SEXP xx, yy, zz, ff = R_NilValue, res = R_NilValue;

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

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    res = ast_discrete(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, test_type, a, debuglevel);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    res = ast_gaustests(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, a, debuglevel, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    res = ast_micg(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
            maxsize, a, debuglevel);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

      res = ast_dperm(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
              maxsize, a, test_type, INT(B),
              IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

      res = ast_gperm(xx, yy, zz, ff, x, y, sx, fixed, nobs, nsx, nf, minsize,
              maxsize, a, test_type, INT(B),
              IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/

  UNPROTECT(3 + (nf > 0));

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  return res;

}/*ALLSUBS_TEST*/

