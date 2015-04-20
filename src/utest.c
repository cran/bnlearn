#include "include/rcore.h"
#include "include/allocations.h"
#include "include/globals.h"
#include "include/tests.h"
#include "include/dataframe.h"
#include "include/covariance.h"

#define DISCRETE_SWAP_X() \
      xdata = VECTOR_ELT(xx, i); \
      xptr = INTEGER(xdata); \
      llx = NLEVELS(xdata);

#define GAUSSIAN_SWAP_X() \
      xptr = REAL(VECTOR_ELT(xx, i)); \
      xm = c_mean(xptr, nobs); \
      xsd = c_var(xptr, xm, nobs);

/* mutual information test, with and without df adjustments. */
static inline double ut_mi(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, int adj) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
void *memp = NULL;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 2 * nobs * c_mi(xptr, llx, yptr, lly, nobs, df, adj);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_MI*/

/* shrinkage mutual information test. */
static inline double ut_shmi(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
void *memp = NULL;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 2 * nobs * c_shmi(xptr, llx, yptr, lly, nobs);
    *df = ((double)(llx - 1) * (double)(lly - 1));
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_SHMI*/

/* Jonckheere-Terpstra test. */
static inline double ut_jt(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
void *memp = NULL;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = c_jt(xptr, llx, yptr, lly, nobs);
    pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_JT*/

/* Pearson's X^2 test, with and without df adjustments. */
static inline double ut_x2(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, int adj) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
void *memp = NULL;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = c_x2(xptr, llx, yptr, lly, nobs, df, adj);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_X2*/

/* linear correlation. */
static inline double ut_cor(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0;
double transform = 0, *xptr = NULL, *yptr = REAL(yy);
double xm = 0, ym = 0, xsd = 0, ysd = 0, statistic = 0;

  /* check that we have enough degrees of freedom. */
  *df = nobs - 2;
  if (*df < 1)
    error("trying to do an independence test with zero degrees of freedom.");

  /* cache mean and variance. */
  ym = c_mean(yptr, nobs);
  ysd = c_var(yptr, ym, nobs);

  for (i = 0; i < ntests; i++) {

    /* no allocations require a vmax{get,set}() call. */
    GAUSSIAN_SWAP_X();
    statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
    transform = fabs(statistic * sqrt(*df) / sqrt(1 - statistic * statistic));
    pvalue[i] = 2 * pt(transform, *df, FALSE, FALSE);

  }/*FOR*/

  return statistic;

}/*UT_COR*/

/* Fisher's Z test. */
static inline double ut_zf(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue) {

int i = 0;
double *xptr = NULL, *yptr = REAL(yy);
double xm = 0, ym = 0, xsd = 0, ysd = 0, statistic = 0;

  /* check that we have enough degrees of freedom. */
  if (nobs - 3 < 1)
    error("trying to do an independence test with zero degrees of freedom.");

  /* cache mean and variance. */
  ym = c_mean(yptr, nobs);
  ysd = c_var(yptr, ym, nobs);

  for (i = 0; i < ntests; i++) {

    /* no allocations require a vmax{get,set}() call. */
    GAUSSIAN_SWAP_X();
    statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
    statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                  sqrt((double)nobs - 3);
    pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

  }/*FOR*/

  return statistic;

}/*UT_ZF*/

/* Gaussian mutual information test. */
static inline double ut_mig(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0;
double *xptr = NULL, *yptr = REAL(yy);
double xm = 0, ym = 0, xsd = 0, ysd = 0, statistic = 0;

  /* cache mean and variance. */
  ym = c_mean(yptr, nobs);
  ysd = c_var(yptr, ym, nobs);
  /* set the degrees of freedom. */
  *df = 1;

  for (i = 0; i < ntests; i++) {

    /* no allocations require a vmax{get,set}() call. */
    GAUSSIAN_SWAP_X();
    statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
    statistic = - nobs * log(1 - statistic * statistic);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  return statistic;

}/*UT_MIG*/

/* shrinkage Gaussian mutual information test. */
static inline double ut_shmig(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0;
double *yptr = REAL(yy), statistic = 0;

  /* set the degrees of freedom. */
  *df = 1;

  for (i = 0; i < ntests; i++) {

    /* no allocations require a vmax{get,set}() call. */
    statistic = c_fast_shcor(REAL(VECTOR_ELT(xx, i)), yptr, &nobs);
    statistic = - nobs * log(1 - statistic * statistic);
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  return statistic;

}/*UT_SHMIG*/

/* conditional linear Gaussian mutual information test. */
static inline double ut_micg(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, xtype = 0, ytype = TYPEOF(yy), llx = 0, lly = 0;
double xm = 0, xsd = 0, ym = 0, ysd = 0, statistic = 0;
void *xptr = NULL, *yptr = NULL, *memp = NULL;
SEXP xdata;

  if (ytype == INTSXP) {

    /* cache the number of levels. */
    lly = NLEVELS(yy);
    yptr = INTEGER(yy);

  }/*THEN*/
  else {

    /* cache mean and variance. */
    yptr = REAL(yy);
    ym = c_mean(yptr, nobs);
    ysd = c_var(yptr, ym, nobs);

  }/*ELSE*/

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xtype = TYPEOF(xdata);

    if ((ytype == INTSXP) && (xtype == INTSXP)) {

      memp = vmaxget();

      /* if both nodes are discrete, the test reverts back to a discrete
       * mutual information test. */
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);
      DISCRETE_SWAP_X();
      statistic = 2 * nobs * c_mi(xptr, llx, yptr, lly, nobs, df, FALSE);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      vmaxset(memp);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      /* if both nodes are continuous, the test reverts back to a Gaussian
       * mutual information test. */
      xptr = REAL(xdata);
      xm = c_mean(xptr, nobs);
      xsd = c_var(xptr, xm, nobs);
      statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
      *df = 1;
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else {

      if (xtype == INTSXP) {

        xptr = INTEGER(xdata);
        llx = NLEVELS(xdata);
        statistic = 2 * nobs * c_micg(yptr, ym, ysd, xptr, llx, nobs);
        *df = llx - 1;
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*THEN*/
      else {

        xptr = REAL(xdata);
        xm = c_mean(xptr, nobs);
        xsd = c_var(xptr, xm, nobs);
        statistic = 2 * nobs * c_micg(xptr, xm, xsd, yptr, lly, nobs);
        *df = lly - 1;
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*ELSE*/

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_MICG*/

/* discrete permutation tests. */
static inline double ut_dperm(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, int type, int B, int a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy);
int llx = 0, lly = NLEVELS(yy);
double statistic = 0;
void *memp = NULL;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    DISCRETE_SWAP_X();
    statistic = 0;
    c_mcarlo(xptr, llx, yptr, lly, nobs, B, &statistic, pvalue + i,
      a, type, df);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_DPERM*/

/* continuous permutation tests. */
static inline double ut_gperm(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, int type, int B, int a) {

int i = 0;
double *yptr = REAL(yy), statistic = 0;
void *memp = NULL;

  for (i = 0; i < ntests; i++) {

    memp = vmaxget();

    statistic = 0;
    c_gauss_mcarlo(REAL(VECTOR_ELT(xx, i)), yptr, nobs, B, pvalue + i,
      a, type, &statistic);

    vmaxset(memp);

  }/*FOR*/

  return statistic;

}/*UT_GPERM*/

/* unconditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
SEXP xx, yy, result;

  /* allocate the return value, which has the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  nobs = length(yy);

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    statistic = ut_mi(xx, yy, nobs, ntests, pvalue, &df,
                  (strcmp(t, "mi-adf") == 0));

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    statistic = ut_shmi(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    statistic = ut_x2(xx, yy, nobs, ntests, pvalue, &df,
                  (strcmp(t, "x2-adf") == 0));

  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    statistic = ut_jt(xx, yy, nobs, ntests, pvalue);

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* linear correlation. */
    statistic = ut_cor(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    statistic = ut_zf(xx, yy, nobs, ntests, pvalue);

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    statistic = ut_mig(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    statistic = ut_shmig(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (strcmp(t, "mi-cg") == 0) {

    /* conditional linear Gaussian mutual information test. */
    statistic = ut_micg(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0;
    double a = (strncmp(t, "smc-", 4) == 0) ? NUM(alpha) : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type))
      statistic = ut_dperm(xx, yy, nobs, ntests, pvalue, &df, type, INT(B), a);
    else
      statistic = ut_gperm(xx, yy, nobs, ntests, pvalue, type, INT(B), a);

  }/*THEN*/

  UNPROTECT(3);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*UTEST*/

