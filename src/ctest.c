#include "include/rcore.h"
#include "include/sets.h"
#include "include/allocations.h"
#include "include/dataframe.h"
#include "include/tests.h"
#include "include/covariance.h"
#include "include/globals.h"

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
    column = (double **) alloc1dpointer(ncols); \
    column[1] = REAL(yy); \
    for (i = 0; i < nsx; i++) \
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

#define GAUSSIAN_CACHE() \
    /* allocate the covariance matrix. */ \
    cov = alloc1dreal(ncols * ncols); \
    basecov = alloc1dreal(ncols * ncols); \
    /* allocate the matrices needed for the SVD decomposition. */ \
    u = alloc1dreal(ncols * ncols); \
    d = alloc1dreal(ncols); \
    vt = alloc1dreal(ncols * ncols); \
    GAUSSIAN_COLUMN_CACHE();

#define GAUSSIAN_PCOR_CACHE() \
        /* extract and de-reference the i-th variable. */ \
        column[0] = REAL(VECTOR_ELT(xx, i)); \
        /* update the corresponding mean in the cache. */ \
        c_update_meanvec(column, mean, 0, nobs); \
        /* update the covariance matrix. */ \
        memcpy(cov, basecov, ncols * ncols * sizeof(double)); \
        c_update_covmat(column, mean, 0, ncols, nobs, cov); \
        /* compute the partial correlation and the test statistic. */ \
        statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);

#define GAUSSIAN_PCOR_NOCACHE() \
      /* extract and de-reference the i-th variable. */ \
      column[0] = REAL(VECTOR_ELT(xx, 0)); \
      /* allocate and compute mean values and the covariance matrix. */ \
      mean = alloc1dreal(ncols); \
      c_meanvec(column, mean, nobs, ncols, 0); \
      c_covmat(column, mean, ncols, nobs, cov, 0); \
      /* compute the partial correlation and the test statistic. */ \
      statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);

/* mutual information test, with and without df adjustments. */
static inline double ct_mi(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, int adj) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
void *memp = NULL;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 2 * nobs * c_cmi(xptr, llx, yptr, lly, zptr, llz,
                               nobs, df, adj);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_MI*/

/* shrinkage mutual information test. */
static inline double ct_mish(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
void *memp = NULL;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 2 * nobs * c_shcmi(xptr, llx, yptr, lly, zptr, llz,
                               nobs, df);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_MISH*/

/* Pearson's X^2 test, with and without df adjustments. */
static inline double ct_x2(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, int adj) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
void *memp = NULL;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = c_cx2(xptr, llx, yptr, lly, zptr, llz, nobs, df, adj);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_X2*/

/* Jonckheere-Terpstra test. */
static inline double ct_jt(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
void *memp = NULL;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);
    pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_JT*/

/* Pearson's linear correlation test. */
static inline double ct_cor(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double transform = 0, **column = NULL, *mean = NULL, statistic = 0;
double *u = NULL, *d = NULL, *vt = NULL, *cov = NULL, *basecov = 0;

  /* check that we have enough degrees of freedom. */
  *df = nobs - ncols;

  if (*df < 1)
    error("trying to do a conditional independence test with zero degrees of freedom.");

  GAUSSIAN_CACHE();

  if (ntests > 1) {

    /* allocate and compute mean values and the covariance matrix. */
    mean = alloc1dreal(ncols);
    c_meanvec(column, mean, nobs, ncols, 1);
    c_covmat(column, mean, ncols, nobs, cov, 1);
    memcpy(basecov, cov, ncols * ncols * sizeof(double));

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */

      GAUSSIAN_PCOR_CACHE();
      transform = fabs(statistic * sqrt(*df) / sqrt(1 - statistic * statistic));
      pvalue[i] = 2 * pt(transform, *df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else {

    GAUSSIAN_PCOR_NOCACHE();
    transform = fabs(statistic * sqrt(*df) / sqrt(1 - statistic * statistic));
    pvalue[0] = 2 * pt(transform, *df, FALSE, FALSE);

  }/*ELSE*/

  return statistic;

}/*CT_COR*/

/* Fisher's Z test. */
static inline double ct_zf(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double **column = NULL, *mean = NULL, *cov = NULL, *basecov = NULL;
double *u = NULL, *d = NULL, *vt = NULL, statistic = 0;

  /* check that we have enough degrees of freedom. */
  if (nobs - 1 - ncols < 1)
    error("trying to do an independence test with zero degrees of freedom.");

  GAUSSIAN_CACHE();

  if (ntests > 1) {

    /* allocate and compute mean values and the covariance matrix. */
    mean = alloc1dreal(ncols);
    c_meanvec(column, mean, nobs, ncols, 1);
    c_covmat(column, mean, ncols, nobs, cov, 1);
    memcpy(basecov, cov, ncols * ncols * sizeof(double));

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */

      GAUSSIAN_PCOR_CACHE();
      statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                    sqrt((double)nobs - 1 - ncols);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else {

    GAUSSIAN_PCOR_NOCACHE();
    statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                  sqrt((double)nobs - 1 - ncols);
    pvalue[0] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

  }/*ELSE*/

  return statistic;

}/*CT_ZF*/

/* Gaussian mutual information test. */
static inline double ct_mig(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double **column = NULL, *mean = NULL, *cov = NULL, *basecov = NULL;
double *u = NULL, *d = NULL, *vt = NULL, statistic = 0;

  *df = 1;

  GAUSSIAN_CACHE();

  if (ntests > 1) {

    /* allocate and compute mean values and the covariance matrix. */
    mean = alloc1dreal(ncols);
    c_meanvec(column, mean, nobs, ncols, 1);
    c_covmat(column, mean, ncols, nobs, cov, 1);
    memcpy(basecov, cov, ncols * ncols * sizeof(double));

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */

      GAUSSIAN_PCOR_CACHE();
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else {

    GAUSSIAN_PCOR_NOCACHE();
    statistic = - nobs * log(1 - statistic * statistic);
    pvalue[0] = pchisq(statistic, *df, FALSE, FALSE);

  }/*ELSE*/

  return statistic;

}/*CT_MIG*/

/* shrinkage Gaussian mutual information test. */
static inline double ct_migsh(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double **column = NULL, *mean = NULL, *cov = NULL;
double *u = NULL, *d = NULL, *vt = NULL, statistic = 0;

  *df = 1;

  /* allocate the covariance matrix. */
  cov = alloc1dreal(ncols * ncols);

  /* allocate the matrices needed for the SVD decomposition. */
  u = alloc1dreal(ncols * ncols);
  d = alloc1dreal(ncols);
  vt = alloc1dreal(ncols * ncols);

  GAUSSIAN_COLUMN_CACHE();

  /* allocate and compute the mean values. */
  mean = alloc1dreal(ncols);
  c_meanvec(column, mean, nobs, ncols, 1);

  for (i = 0; i < ntests; i++) {

    /* no allocations require a vmax{get,set}() call. */

    /* extract and de-reference the i-th variable. */
    column[0] = REAL(VECTOR_ELT(xx, i));
    /* update the corresponding mean in the cache. */
    c_update_meanvec(column, mean, 0, nobs);
    /* compute the covariance matrix. */
    memset(cov, '\0', ncols * ncols * sizeof(double));
    c_cov_lambda(column, mean, ncols, nobs, cov);
    /* compute the partial correlation and the test statistic. */
    statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
    statistic = - nobs * log(1 - statistic * statistic);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  return statistic;

}/*CT_MIGSH*/


/* conditional linear Gaussian mutual information test. */
static inline double ct_micg(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df) {

int xtype = 0, ytype = TYPEOF(yy), *nlvls = NULL, llx = 0, lly = 0, llz = 0;
int ndp = 0, ngp = 0, nsx = length(zz), **dp = NULL, *dlvls = NULL, j = 0, k = 0;
int i = 0, *zptr = 0;
void *xptr = NULL, *yptr = NULL, **columns = NULL;
double **gp = NULL;
double statistic = 0;
SEXP xdata, temp;

  if (ytype == INTSXP) {

    /* cache the number of levels. */
    lly = NLEVELS(yy);
    yptr = INTEGER(yy);

  }/*THEN*/
  else {

    yptr = REAL(yy);

  }/*ELSE*/

  /* extract the conditioning variables and cache their types. */
  columns = Calloc(nsx, void *);
  nlvls = Calloc(nsx, int);

  for (i = 0; i < nsx; i++) {

    temp = VECTOR_ELT(zz, i);

    if (TYPEOF(temp) == INTSXP) {

      ndp++;
      columns[i] = INTEGER(temp);
      nlvls[i] = NLEVELS(temp);

    }/*THEN*/
    else {

      ngp++;
      columns[i] = REAL(temp);

    }/*ELSE*/

  }/*FOR*/

  dp = Calloc(ndp + 1, int *);
  gp = Calloc(ngp + 1, double *);
  dlvls = Calloc(ndp + 1, int);
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

    zptr = Calloc(nobs, int);
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
        statistic = 2 * nobs * c_cmi(xptr, llx, yptr, lly, zptr, llz,
                                 nobs, df, FALSE);

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

  Free(columns);
  Free(nlvls);
  Free(dlvls);
  Free(zptr);
  Free(dp);
  Free(gp);

  return statistic;

}/*CT_MICG*/

/* discrete permutation tests. */
static inline double ct_dperm(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, int type, int B, int a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
int llx = 0, lly = NLEVELS(yy), llz = 0;
double statistic = 0;
void *memp = NULL;
SEXP xdata, config;

  DISCRETE_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 0;
    c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, B, &statistic,
      pvalue + i, a, type, df);

    vmaxset(memp);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_DPERM*/

/* continuous permutation tests. */
static inline double ct_gperm(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, int type, int B, int a) {

int i = 0, nsx = length(zz), ncols = nsx + 2;
double **column = NULL, *yptr = REAL(yy), statistic = 0;
void *memp = NULL;

  GAUSSIAN_COLUMN_CACHE();

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    /* swap the first column and restore the second, which is that undergoing
     * permutation (backward compatibility from set random seed). */
    column[0] = REAL(VECTOR_ELT(xx, i));
    column[1] = yptr;

    statistic = 0;
    c_gauss_cmcarlo(column, ncols, nobs, B, &statistic, pvalue + i,
      a, type);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*CT_GPERM*/

/* conditional independence tests. */
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
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

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    statistic = ct_mi(xx, yy, zz, nobs, ntests, pvalue, &df,
                  (strcmp(t, "mi-adf") == 0));

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    statistic = ct_mish(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    statistic = ct_x2(xx, yy, zz, nobs, ntests, pvalue, &df,
                  (strcmp(t, "x2-adf") == 0));

  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    statistic = ct_jt(xx, yy, zz, nobs, ntests, pvalue);

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* Pearson's linear correlation test. */
    statistic = ct_cor(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    statistic = ct_zf(xx, yy, zz, nobs, ntests, pvalue);

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    statistic = ct_mig(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    statistic = ct_migsh(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "mi-cg") == 0) {

    /* conditional linear Gaussian mutual information test. */
    statistic = ct_micg(xx, yy, zz, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0;
    double a = (strncmp(t, "smc-", 4) == 0) ? NUM(alpha) : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type))
      statistic = ct_dperm(xx, yy, zz, nobs, ntests, pvalue, &df, type, INT(B), a);
    else
      statistic = ct_gperm(xx, yy, zz, nobs, ntests, pvalue, &df, type, INT(B), a);

  }/*THEN*/

  UNPROTECT(4);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*CTEST*/

