#include "include/rcore.h"
#include "include/sets.h"
#include "include/dataframe.h"
#include "include/tests.h"
#include "include/covariance.h"
#include "include/globals.h"
#include "include/blas.h"

#define DISCRETE_SWAP_X() \
      xdata = VECTOR_ELT(xx, i); \
      xptr = INTEGER(xdata); \
      llx = NLEVELS(xdata);

#define DISCRETE_CACHE() \
    PROTECT(config = c_configurations(zz, TRUE, TRUE)); \
    zptr = INTEGER(config); \
    llz = NLEVELS(config);

#define GAUSSIAN_COLUMN_CACHE() \
    /* allocate and initialize an array of pointers for the variables. */ \
    column = Calloc1D(ncols, sizeof(double *)); \
    column[1] = REAL(yy); \
    for (i = 0; i < nsx; i++) \
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

#define GAUSSIAN_CACHE() \
    /* allocate the covariance matrix. */ \
    cov = Calloc1D(ncols * ncols, sizeof(double)); \
    basecov = Calloc1D(ncols * ncols, sizeof(double)); \
    /* allocate the matrices needed for the SVD decomposition. */ \
    c_udvt(&u, &d, &vt, ncols); \
    GAUSSIAN_COLUMN_CACHE();

#define GAUSSIAN_PCOR_CACHE() \
        /* extract and de-reference the i-th variable. */ \
        column[0] = REAL(VECTOR_ELT(xx, i)); \
        /* update the corresponding mean in the cache. */ \
        c_update_meanvec(column, mean, 0, nobs); \
        /* update the covariance matrix. */ \
        memcpy(cov, basecov, ncols * ncols * sizeof(double)); \
        c_update_covmat(column, mean, 0, ncols, nobs, cov); \

#define GAUSSIAN_PCOR_NOCACHE() \
      /* extract and de-reference the i-th variable. */ \
      column[0] = REAL(VECTOR_ELT(xx, 0)); \
      /* allocate and compute mean values and the covariance matrix. */ \
      mean = Calloc1D(ncols, sizeof(double)); \
      c_meanvec(column, mean, nobs, ncols, 0); \
      c_covmat(column, mean, ncols, nobs, cov, 0);

#define COMPUTE_PCOR() \
        /* compute the partial correlation and the test statistic. */ \
        statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);

#define GAUSSIAN_FREE() \
  Free1D(u); \
  Free1D(vt); \
  Free1D(d); \
  Free1D(cov); \
  Free1D(basecov);

/* parametric tests for discrete variables. */
static double ct_discrete(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, nobs, df, test);
      if ((test == MI) || (test == MI_ADF))
        statistic = 2 * nobs * statistic;
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = 2 * nobs * c_shcmi(xptr, llx, yptr, lly, zptr, llz,
                               nobs, df);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_DISCRETE*/

/* parametric tests for Gaussian variables. */
static double ct_gaustests(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double transform = 0, **column = NULL, *mean = NULL, statistic = 0, lambda = 0;
double *u = NULL, *d = NULL, *vt = NULL, *cov = NULL, *basecov = 0;

  /* compute the degrees of freedom for correlation and mutual information. */
  if (test == COR)
    *df = nobs - ncols;
  else if ((test == MI_G) || (test == MI_G_SH))
    *df = 1;

  if (((test == COR) && (*df < 1)) || ((test == ZF) && (nobs - ncols < 2))) {

    /* if there are not enough degrees of freedom, return independence. */
    warning("trying to do a conditional independence test with zero degrees of freedom.");

    *df = 0;
    statistic = 0;
    for (i = 0; i < ntests; i++)
      pvalue[i] = 1;

    return statistic;

  }/*THEN*/

  GAUSSIAN_CACHE();

  if (ntests > 1) {

    /* allocate and compute mean values and the covariance matrix. */
    mean = Calloc1D(ncols, sizeof(double));
    c_meanvec(column, mean, nobs, ncols, 1);
    c_covmat(column, mean, ncols, nobs, cov, 1);
    memcpy(basecov, cov, ncols * ncols * sizeof(double));

    for (i = 0; i < ntests; i++) {

      GAUSSIAN_PCOR_CACHE();

      if (test == COR) {

        COMPUTE_PCOR();
        transform = cor_t_trans(statistic, *df);
        pvalue[i] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

      }/*THEN*/
      else if (test == MI_G) {

        COMPUTE_PCOR();
        statistic = 2 * nobs * cor_mi_trans(statistic);
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*THEN*/
      else if (test == MI_G_SH) {

        lambda = covmat_lambda(column, mean, cov, nobs, ncols);
        covmat_shrink(cov, ncols, lambda);
        COMPUTE_PCOR();
        statistic = 2 * nobs * cor_mi_trans(statistic);
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*THEN*/
      else if (test == ZF) {

        COMPUTE_PCOR();
        statistic = cor_zf_trans(statistic, (double)nobs - ncols);
        pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

      }/*THEN*/

    }/*FOR*/

  }/*THEN*/
  else {

    GAUSSIAN_PCOR_NOCACHE();

    if (test == COR) {

      COMPUTE_PCOR();
      transform = cor_t_trans(statistic, *df);
      pvalue[0] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      COMPUTE_PCOR();
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[0] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(column, mean, cov, nobs, ncols);
      covmat_shrink(cov, ncols, lambda);
      COMPUTE_PCOR();
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[0] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      COMPUTE_PCOR();
      statistic = cor_zf_trans(statistic, (double)nobs - ncols);
      pvalue[0] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*ELSE*/

  GAUSSIAN_FREE();

  Free1D(mean);
  Free1D(column);

  return statistic;

}/*CT_GAUSTESTS*/

/* conditional linear Gaussian mutual information test. */
static double ct_micg(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int xtype = 0, ytype = TYPEOF(yy), *nlvls = NULL, llx = 0, lly = 0, llz = 0;
int ndp = 0, ngp = 0, nsx = length(zz), **dp = NULL, *dlvls = NULL, j = 0, k = 0;
int i = 0, *zptr = 0;
void *xptr = NULL, *yptr = NULL, **columns = NULL;
double **gp = NULL;
double statistic = 0;
SEXP xdata;

  if (ytype == INTSXP) {

    /* cache the number of levels. */
    lly = NLEVELS(yy);
    yptr = INTEGER(yy);

  }/*THEN*/
  else {

    yptr = REAL(yy);

  }/*ELSE*/

  /* extract the conditioning variables and cache their types. */
  columns = Calloc1D(nsx, sizeof(void *));
  nlvls = Calloc1D(nsx, sizeof(int));
  df2micg(zz, columns, nlvls, &ndp, &ngp);

  dp = Calloc1D(ndp + 1, sizeof(int *));
  gp = Calloc1D(ngp + 1, sizeof(double *));
  dlvls = Calloc1D(ndp + 1, sizeof(int));
  for (i = 0, j = 0, k = 0; i < nsx; i++)
    if (nlvls[i] > 0) {

      dlvls[1 + j] = nlvls[i];
      dp[1 + j++] = columns[i];

    }/*THEN*/
    else {

      gp[1 + k++] = columns[i];

    }/*ELSE*/

  /* allocate vector for the configurations of the discrete parents; or, if
   * there no discrete parents, for the means of the continuous parents. */
  if (ndp > 0) {

    zptr = Calloc1D(nobs, sizeof(int));
    c_fast_config(dp + 1, nobs, ndp, dlvls + 1, zptr, &llz, 1);

  }/*THEN*/

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xtype = TYPEOF(xdata);

    if (xtype == INTSXP) {

      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

    }/*THEN*/
    else {

      xptr = REAL(xdata);

    }/*ELSE*/

    if ((ytype == INTSXP) && (xtype == INTSXP)) {

      if (ngp > 0) {

        /* need to reverse conditioning to actually compute the test. */
        statistic = 2 * nobs * nobs *
                      c_cmicg_unroll(xptr, llx, yptr, lly, zptr, llz,
                                 gp + 1, ngp, df, nobs);

      }/*THEN*/
      else {

        /* the test reverts back to a discrete mutual information test. */
        statistic = 2 * nobs * c_cchisqtest(xptr, llx, yptr, lly, zptr, llz,
                                 nobs, df, MI);

      }/*ELSE*/

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      gp[0] = xptr;
      statistic = 2 * nobs * c_cmicg(yptr, gp, ngp + 1, NULL, 0, zptr, llz,
                               dlvls, nobs);
      /* one regression coefficient for each conditioning level is added;
       * if all conditioning variables are continuous that's just one global
       * regression coefficient. */
      *df = (llz == 0) ? 1 : llz;

    }/*THEN*/
    else if ((ytype == INTSXP) && (xtype == REALSXP)) {

      dp[0] = yptr;
      dlvls[0] = lly;
      statistic = 2 * nobs * c_cmicg(xptr, gp + 1, ngp, dp, ndp + 1, zptr,
                               llz, dlvls, nobs);

      /* for each additional configuration of the discrete conditioning
       * variables plus the discrete yptr, one whole set of regression
       * coefficients (plus the intercept) is added. */
      *df = (lly - 1) * ((llz == 0) ? 1 : llz)  * (ngp + 1);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == INTSXP)) {

      dp[0] = xptr;
      dlvls[0] = llx;
      statistic = 2 * nobs * c_cmicg(yptr, gp + 1, ngp, dp, ndp + 1, zptr,
                               llz, dlvls, nobs);
      /* same as above, with xptr and yptr swapped. */
      *df = (llx - 1) * ((llz == 0) ? 1 : llz) * (ngp + 1);

    }/*ELSE*/

    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  Free1D(columns);
  Free1D(nlvls);
  Free1D(dlvls);
  Free1D(zptr);
  Free1D(dp);
  Free1D(gp);

  return statistic;

}/*CT_MICG*/

/* discrete permutation tests. */
static double ct_dperm(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e type, int B, double a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
int llx = 0, lly = NLEVELS(yy), llz = 0;
double statistic = 0;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();
    statistic = 0;
    c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, B, &statistic,
      pvalue + i, a, type, df);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_DPERM*/

/* continuous permutation tests. */
static double ct_gperm(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e type, int B, double a) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double **column = NULL, *yptr = REAL(yy), statistic = 0;

  GAUSSIAN_COLUMN_CACHE();

  for (i = 0; i < ntests; i++) {

    /* swap the first column and restore the second, which is that undergoing
     * permutation (backward compatibility from set random seed). */
    column[0] = REAL(VECTOR_ELT(xx, i));
    column[1] = yptr;

    statistic = 0;
    c_gauss_cmcarlo(column, ncols, nobs, B, &statistic, pvalue + i,
      a, type);

  }/*FOR*/

  Free1D(column);

  return statistic;

}/*CT_GPERM*/

/* conditional independence tests. */
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
SEXP xx, yy, zz, result;

  /* allocate the return value, which the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, FALSE));
  nobs = length(yy);

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    statistic = ct_discrete(xx, yy, zz, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    statistic = ct_gaustests(xx, yy, zz, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    statistic = ct_micg(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    statistic = ct_dperm(xx, yy, zz, nobs, ntests, pvalue, &df, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    statistic = ct_gperm(xx, yy, zz, nobs, ntests, pvalue, &df, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/

  UNPROTECT(4);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statstic '%s'.", t);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*CTEST*/

