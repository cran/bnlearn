#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/dataframe.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"

#define SKIP_FIXED() \
    /* nothing to do if there are no more valid nodes.*/ \
    if (valid - 1 == 0) break; \
    /* do not test fixed variables. */ \
    if (ff[i] > 0) continue;

#define DISCRETE_CACHE() \
  /* allocate and initialize an array of pointers for the variables. */ \
  column = (int **) Calloc1D(nz, sizeof(int *)); \
  for (i = 0; i < nz; i++) \
    column[i] = INTEGER(VECTOR_ELT(zz, i)); \
  /* allocate and compute the number of levels. */ \
  nlvls = Calloc1D(nz, sizeof(int)); \
  for (i = 0; i < nz; i++) \
    nlvls[i] = NLEVELS(VECTOR_ELT(zz, i)); \
  /* allocate the parents' configurations. */ \
  zptr = Calloc1D(nobs, sizeof(int));

#define FREE_DISCRETE_CACHE() \
  Free1D(column); \
  Free1D(nlvls); \
  Free1D(zptr);

#define GAUSSIAN_CACHE() \
  /* allocate and initialize an array of pointers for the variables. */ \
  column = (double **) Calloc1D(nz, sizeof(double *)); \
  for (i = 0; i < nz; i++) \
    column[i] = REAL(VECTOR_ELT(zz, i)); \
  /* allocate and compute mean values and the covariance matrix. */ \
  mean = Calloc1D(nz, sizeof(double)); \
  c_meanvec(column, mean, nobs, nz, 0);

#define FREE_GAUSSIAN_CACHE() \
  Free1D(column); \
  Free1D(mean);

#define GAUSSIAN_ALLOC() \
  cov = Calloc1D((nz + 1) * (nz + 1), sizeof(double)); \
  c_udvt(&u, &d, &vt, nz + 1);

#define FREE_GAUSSIAN_ALLOC() \
  Free1D(cov);

#define ALLOC_DISCRETE_SUBSET() \
 /* allocate the subset. */ \
 subset = (int **) Calloc1D(nz - 1, sizeof(int *)); \
 sublvls = Calloc1D(nz - 1, sizeof(int));

#define FREE_DISCRETE_SUBSET() \
  Free1D(subset); \
  Free1D(sublvls);

#define ALLOC_GAUSSIAN_SUBSET() \
  /* allocate the subset. */ \
  subset = (double **) Calloc1D(nz + 1, sizeof(double *)); \
  subset[0] = xptr; \
  submean = Calloc1D(nz + 1, sizeof(double)); \
  submean[0] = c_mean(xptr, nobs);

#define FREE_GAUSSIAN_SUBSET() \
  Free1D(subset); \
  Free1D(submean);

#define PREPARE_DISCRETE_SUBSET() \
    /* copy all valid variables in the subset but the current one. */ \
    for (j = 0, k = 0; j < nz; j++) { \
      if (column[j] == NULL) \
        continue; \
      if (j == i) { \
        yptr = column[j]; \
        lly = nlvls[j]; \
      }/*THEN*/ \
      else { \
        subset[k] = column[j]; \
        sublvls[k++] = nlvls[j]; \
      }/*ELSE*/ \
    }/*FOR*/ \
    /* construct the parents' configurations. */ \
    c_fast_config(subset, nobs, valid - 1, sublvls, zptr, &llz, 1);

#define PREPARE_GAUSSIAN_SUBSET() \
    /* copy all valid variables in the subset but the current one. */ \
    for (j = 0, k = 2; j < nz; j++) { \
      if (column[j] == NULL) \
        continue; \
      if (j == i) { \
        subset[1] = column[j]; \
        submean[1] = mean[j]; \
      }/*THEN*/ \
      else { \
        subset[k] = column[j]; \
        submean[k++] = mean[j]; \
      }/*ELSE*/ \
    }/*FOR*/

#define FLAG_AND_DEBUG() \
  if (pvalue[cur - 1] > alpha) { \
    if (debuglevel > 0) { \
      Rprintf("    > node %s is independent from %s given ", \
        CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(z, i))); \
      for (j = 0; j < nz; j++) \
        if ((j != i) && (column[j] != NULL)) \
          Rprintf("%s ", CHAR(STRING_ELT(z, j))); \
      Rprintf("(p-value: %g).\n", pvalue[cur - 1]); \
    }/*THEN*/ \
    /* remove from the set of valid candidates for the conditioning set. */ \
    valid--; \
    column[cur - 1] = NULL; \
  }/*THEN*/ \
  else { \
    if (debuglevel > 0) { \
      Rprintf("    > node %s is dependent on %s given ", \
        CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(z, i))); \
      for (j = 0; j < nz; j++) \
        if ((j != i) && (column[j] != NULL)) \
          Rprintf("%s ", CHAR(STRING_ELT(z, j))); \
      Rprintf("(p-value: %g).\n", pvalue[cur - 1]); \
    }/*THEN*/ \
  }/*ELSE*/

#define GAUSSIAN_FREE() \
  Free1D(u); \
  Free1D(vt); \
  Free1D(d);

/* parametric tests for discrete variables. */
static void rrd_discrete(SEXP xx, SEXP zz, int *ff, SEXP x, SEXP z, int nobs,
    int nz, test_e test, double *pvalue, double alpha, int debuglevel) {

int i = 0, j = 0, k = 0, cur = 0, llx = NLEVELS(xx), lly = 0, llz = 0;
int *xptr = INTEGER(xx), *yptr = NULL, *zptr = NULL, valid = nz;
int **column = NULL, **subset = NULL, *sublvls = NULL, *nlvls = NULL;
double statistic = 0, df = 0;

  DISCRETE_CACHE();
  ALLOC_DISCRETE_SUBSET();

  for (i = 0; i < nz; i++) {

    SKIP_FIXED();
    PREPARE_DISCRETE_SUBSET();

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, nobs, &df, test);
      if ((test == MI) || (test == MI_ADF))
        statistic = 2 * nobs * statistic;
      pvalue[cur++] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, nobs, &df);
      pvalue[cur++] = pchisq(2 * nobs * statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
      pvalue[cur++] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    FLAG_AND_DEBUG();

  }/*FOR*/

  FREE_DISCRETE_CACHE();
  FREE_DISCRETE_SUBSET();

}/*RRD_DISCRETE*/

/* parametric tests for discrete variables. */
static void rrd_gaustests(SEXP xx, SEXP zz, int *ff, SEXP x, SEXP z, int nobs,
    int nz, test_e test, double *pvalue, double alpha, int debuglevel) {

int i = 0, j = 0, k = 0, cur = 0, valid = nz, df = 0;
double *xptr = REAL(xx);
double **column = NULL, **subset = NULL, *submean = NULL, *mean = NULL;
double *cov = NULL, *u = NULL, *d = NULL, *vt = NULL, lambda = 0;
double statistic = 0;

  GAUSSIAN_CACHE();
  GAUSSIAN_ALLOC();
  ALLOC_GAUSSIAN_SUBSET();

  for (i = 0; i < nz; i++) {

    SKIP_FIXED();

    /* compute the degrees of freedom for correlation and mutual information. */
    if (test == COR)
      df = nobs - (valid + 1);
    else if ((test == MI_G) || (test == MI_G_SH))
      df = 1;

    if (((test == COR) && (df < 1)) ||
        ((test == ZF) && (nobs - (valid + 1) < 2))) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      pvalue[cur++] = 1;
      FLAG_AND_DEBUG();
      continue;

    }/*THEN*/

    PREPARE_GAUSSIAN_SUBSET();

    /* compute the covariance matrix. */
    c_covmat(subset, submean, valid + 1, nobs, cov, 0);

    if (test == COR) {

      statistic = c_fast_pcor(cov, u, d, vt, valid + 1, TRUE);
      statistic = cor_t_trans(statistic, (double)df);
      pvalue[cur++] = 2 * pt(fabs(statistic), df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = c_fast_pcor(cov, u, d, vt, valid + 1, TRUE);
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[cur++] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(subset, submean, cov, nobs, valid + 1);
      covmat_shrink(cov, valid + 1, lambda);
      statistic = c_fast_pcor(cov, u, d, vt, valid + 1, TRUE);
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[cur++] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = c_fast_pcor(cov, u, d, vt, valid + 1, TRUE);
      statistic = cor_zf_trans(statistic, (double)nobs - (valid + 1));
      pvalue[cur++] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    FLAG_AND_DEBUG();

  }/*FOR*/

  GAUSSIAN_FREE();
  FREE_GAUSSIAN_CACHE();
  FREE_GAUSSIAN_ALLOC();
  FREE_GAUSSIAN_SUBSET();

}/*RRD_GAUSTESTS*/

/* conditional linear Gaussian test. */
static void rrd_micg(SEXP xx, SEXP zz, int *ff, SEXP x, SEXP z, int nobs,
    int nz, test_e test, double *pvalue, double alpha, int debuglevel) {

int i = 0, j = 0, k = 0, l = 0, ndp = 0, ngp = 0;
int cur = 0, valid = nz, xtype = TYPEOF(xx), llx = 0, lly = 0, llz = 0;
int *nlvls = NULL, *zptr = NULL, **dp = NULL, *dlvls = NULL;
double statistic = 0, df = 0, **gp = NULL;
void *xptr = NULL, *yptr = NULL, **column = NULL;

  /* scan the parents for type and number of levels. */
  column = Calloc1D(nz, sizeof(void *));
  nlvls = Calloc1D(nz, sizeof(int));
  df2micg(zz, column, nlvls, &ndp, &ngp);

  /* dereference xx, whatever the type. */
  if (xtype == INTSXP) {

    xptr = INTEGER(xx);
    llx = NLEVELS(xx);

  }/*THEN*/
  else {

    xptr = REAL(xx);

  }/*ELSE*/

  for (i = 0; i < nz; i++) {

    SKIP_FIXED();

    /* count how many discrete and continuous parents are in the subset. */
    for (j = 0, ndp = 0, ngp = 0; j < nz; j++) {

      if ((column[j] == NULL) || (i == j))
        continue;

      if (nlvls[j] > 0)
        ndp++;
      else
        ngp++;

    }/*FOR*/

    /* split parents by type. */
    gp = Calloc1D(ngp + 1, sizeof(double *));
    dp = Calloc1D(ndp + 1, sizeof(int *));
    dlvls = Calloc1D(ndp + 1, sizeof(int));
    for (j = 0, k = 0, l = 0; j < nz; j++) {

      if ((column[j] == NULL) || (i == j))
        continue;

      if (nlvls[j] > 0) {

        dp[1 + k] = column[j];
        dlvls[1 + k++] = nlvls[j];

      }/*THEN*/
      else {

        gp[1 + l++] = column[j];

      }/*ELSE*/

    }/*FOR*/

    /* extract the target variable. */
    if (nlvls[i] > 0) {

      yptr = column[i];
      lly = nlvls[i];

    }/*THEN*/
    else {

      yptr = column[i];

    }/*ELSE*/

    /* if there are discrete conditioning variables, compute their
     * configurations. */
    if (ndp  > 0) {

      zptr = Calloc1D(nobs, sizeof(int));
      c_fast_config(dp + 1, nobs, ndp, dlvls + 1, zptr, &llz, 1);

    }/*THEN*/
    else {

      zptr = NULL;
      llz = 0;

    }/*ELSE*/

    if ((xtype == INTSXP) && (nlvls[i] > 0)) {

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
    else if ((xtype == REALSXP) && (nlvls[i] == 0)) {

      gp[0] = xptr;
      statistic = 2 * nobs * c_cmicg(yptr, gp, ngp + 1, NULL, 0, zptr, llz,
                             dlvls, nobs);

      /* one regression coefficient for each conditioning level is added;
       * if all conditioning variables are continuous that's just one global
       * regression coefficient. */
      df = (llz == 0) ? 1 : llz;

    }/*THEN*/
    else if ((xtype == REALSXP) && (nlvls[i] > 0)) {

      dp[0] = yptr;
      dlvls[0] = lly;
      statistic = 2 * nobs * c_cmicg(xptr, gp + 1, ngp, dp, ndp + 1, zptr,
                               llz, dlvls, nobs);

      /* for each additional configuration of the discrete conditioning
       * variables plus the discrete yptr, one whole set of regression
       * coefficients (plus the intercept) is added. */
      df = (lly - 1) * ((llz == 0) ? 1 : llz)  * (ngp + 1);

    }/*THEN*/
    else if ((xtype == INTSXP) && (nlvls[i] == 0)) {

        dp[0] = xptr;
        dlvls[0] = llx;
        statistic = 2 * nobs * c_cmicg(yptr, gp + 1, ngp, dp, ndp + 1, zptr,
                                 llz, dlvls, nobs);

        /* for each additional configuration of the discrete conditioning
         * variables plus the discrete yptr, one whole set of regression
         * coefficients (plus the intercept) is added. */
        df = (llx - 1) * ((llz == 0) ? 1 : llz)  * (ngp + 1);

    }/*THEN*/

    pvalue[cur++] = pchisq(statistic, df, FALSE, FALSE);

    Free1D(gp);
    Free1D(dp);
    Free1D(dlvls);
    Free1D(zptr);

    FLAG_AND_DEBUG();

  }/*FOR*/

  Free1D(column);
  Free1D(nlvls);

}/*RRD_MICG*/

/* discrete permutation tests. */
static void rrd_dperm(SEXP xx, SEXP zz, int *ff, SEXP x, SEXP z, int nobs,
    int nz, test_e test, double *pvalue, double alpha, int nperms,
    double threshold, int debuglevel) {

int i = 0, j = 0, k = 0, cur = 0, llx = NLEVELS(xx), lly = 0, llz = 0;
int *xptr = INTEGER(xx), *yptr = NULL, *zptr = NULL, valid = nz;
int **column = NULL, **subset = NULL, *sublvls = NULL, *nlvls = NULL;
double statistic = 0, df = 0;

  DISCRETE_CACHE();
  ALLOC_DISCRETE_SUBSET();

  for (i = 0; i < nz; i++) {

    SKIP_FIXED();
    PREPARE_DISCRETE_SUBSET();

    c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, nperms, &statistic,
      pvalue + cur, threshold, test, &df);
    cur++;

    FLAG_AND_DEBUG();

  }/*FOR*/

  FREE_DISCRETE_CACHE();
  FREE_DISCRETE_SUBSET();

}/*RRD_DPERM*/

/* continuous permutation tests. */
static void rrd_gperm(SEXP xx, SEXP zz, int *ff, SEXP x, SEXP z, int nobs,
    int nz, test_e test, double *pvalue, double alpha, int nperms,
    double threshold, int debuglevel) {

int i = 0, j = 0, k = 0, cur = 0, valid = nz;
double *xptr = REAL(xx);
double **column = NULL, **subset = NULL, *submean = NULL, *mean = NULL;
double statistic = 0;

  GAUSSIAN_CACHE();
  ALLOC_GAUSSIAN_SUBSET();

  for (i = 0; i < nz; i++) {

    SKIP_FIXED();
    PREPARE_GAUSSIAN_SUBSET();

    c_gauss_cmcarlo(subset, valid + 1, nobs, nperms, &statistic, pvalue + cur,
      threshold, test);
    cur++;

    FLAG_AND_DEBUG();

  }/*FOR*/

  FREE_GAUSSIAN_CACHE();
  FREE_GAUSSIAN_SUBSET();

}/*RRD_GPERM*/

SEXP roundrobin_test(SEXP x, SEXP z, SEXP fixed, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP debug) {

int nz = length(z), nf = length(fixed), nobs = 0;
int *ff = NULL, debuglevel = isTRUE(debug);
double *pvalue = NULL, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
SEXP xx, zz, try, result;

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, z, FALSE, FALSE));
  nobs = length(xx);
  /* match fixed variables. */
  PROTECT(try = match(fixed, z, 0));
  ff = INTEGER(try);
  /* allocate the return value. */
  PROTECT(result = allocVector(REALSXP, nz - nf));
  setAttrib(result, R_NamesSymbol, string_setdiff(z, fixed));
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', (nz - nf) * sizeof(double));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    rrd_discrete(xx, zz, ff, x, z, nobs, nz, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    rrd_gaustests(xx, zz, ff, x, z, nobs, nz, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    rrd_micg(xx, zz, ff, x, z, nobs, nz, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    rrd_dperm(xx, zz, ff, x, z, nobs, nz, test_type, pvalue, a, INT(B),
      IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    rrd_gperm(xx, zz, ff, x, z, nobs, nz, test_type, pvalue, a, INT(B),
      IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/

  UNPROTECT(4);

  return result;

}/*ROUNDROBIN_TEST*/
