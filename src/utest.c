#include "include/rcore.h"
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
      xsd = c_sse(xptr, xm, nobs);

/* parametric tests for discrete variables. */
static double ut_discrete(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_chisqtest(xptr, llx, yptr, lly, nobs, df, test);
      if ((test == MI) || (test == MI_ADF))
        statistic = 2 * nobs * statistic;
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = 2 * nobs * c_shmi(xptr, llx, yptr, lly, nobs);
      *df = ((double)(llx - 1) * (double)(lly - 1));
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_jt(xptr, llx, yptr, lly, nobs);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_DISCRETE*/

/* parametric tests for Gaussian variables. */
static double ut_gaustests(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0;
double transform = 0, *xptr = NULL, *yptr = REAL(yy);
double xm = 0, ym = 0, xsd = 0, ysd = 0, statistic = 0;

  /* compute the degrees of freedom for correlation and mutual information. */
  if (test == COR)
    *df = nobs - 2;
  else if ((test == MI_G) || (test == MI_G_SH))
    *df = 1;

  if (((test == COR) && (*df < 1)) || ((test == ZF) && (nobs - 2 < 2))) {

    /* if there are not enough degrees of freedom, return independence. */
    warning("trying to do an independence test with zero degrees of freedom.");

    *df = 0;
    statistic = 0;
    for (i = 0; i < ntests; i++)
      pvalue[i] = 1;

     return statistic;

  }/*THEN*/

  /* cache mean and variance. */
  ym = c_mean(yptr, nobs);
  ysd = c_sse(yptr, ym, nobs);

  for (i = 0; i < ntests; i++) {

    GAUSSIAN_SWAP_X();
    statistic = c_fast_cor(xptr, yptr, nobs, xm, ym, xsd, ysd);

    if (test == COR) {

      transform = cor_t_trans(statistic, *df);
      pvalue[i] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      statistic *= 1 - cor_lambda(xptr, yptr, nobs, xm, ym, xsd, ysd, statistic);
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = cor_zf_trans(statistic, (double)nobs - 2);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_GAUSTESTS*/

/* conditional linear Gaussian mutual information test. */
static double ut_micg(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
  double *df) {

int i = 0, xtype = 0, ytype = TYPEOF(yy), llx = 0, lly = 0;
double xm = 0, xsd = 0, ym = 0, ysd = 0, statistic = 0;
void *xptr = NULL, *yptr = NULL;
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
    ysd = c_sse(yptr, ym, nobs);

  }/*ELSE*/

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xtype = TYPEOF(xdata);

    if ((ytype == INTSXP) && (xtype == INTSXP)) {

      /* if both nodes are discrete, the test reverts back to a discrete
       * mutual information test. */
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);
      DISCRETE_SWAP_X();
      statistic = 2 * nobs * c_chisqtest(xptr, llx, yptr, lly, nobs, df, MI);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      /* if both nodes are continuous, the test reverts back to a Gaussian
       * mutual information test. */
      xptr = REAL(xdata);
      xm = c_mean(xptr, nobs);
      xsd = c_sse(xptr, xm, nobs);
      statistic = c_fast_cor(xptr, yptr, nobs, xm, ym, xsd, ysd);
      *df = 1;
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else {

      if (xtype == INTSXP) {

        xptr = INTEGER(xdata);
        llx = NLEVELS(xdata);
        ysd = sqrt(ysd / (nobs - 1));
        statistic = 2 * nobs * c_micg(yptr, ym, ysd, xptr, llx, nobs);
        *df = llx - 1;
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*THEN*/
      else {

        xptr = REAL(xdata);
        xm = c_mean(xptr, nobs);
        xsd = sqrt(c_sse(xptr, xm, nobs) / (nobs - 1));
        statistic = 2 * nobs * c_micg(xptr, xm, xsd, yptr, lly, nobs);
        *df = lly - 1;
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*ELSE*/

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_MICG*/

/* discrete permutation tests. */
static double ut_dperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df, test_e type, int B, double a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy);
int llx = 0, lly = NLEVELS(yy);
double statistic = 0;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();
    statistic = 0;
    c_mcarlo(xptr, llx, yptr, lly, nobs, B, &statistic, pvalue + i,
      a, type, df);

  }/*FOR*/

  return statistic;

}/*UT_DPERM*/

/* continuous permutation tests. */
static double ut_gperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    test_e type, int B, double a) {

int i = 0;
double *yptr = REAL(yy), statistic = 0;

  for (i = 0; i < ntests; i++) {

    statistic = 0;
    c_gauss_mcarlo(REAL(VECTOR_ELT(xx, i)), yptr, nobs, B, pvalue + i,
      a, type, &statistic);

  }/*FOR*/

  return statistic;

}/*UT_GPERM*/

/* unconditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
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

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    statistic = ut_discrete(xx, yy, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    statistic = ut_gaustests(xx, yy, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    statistic = ut_micg(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    statistic = ut_dperm(xx, yy, nobs, ntests, pvalue, &df, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    statistic = ut_gperm(xx, yy, nobs, ntests, pvalue, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/

  UNPROTECT(3);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*UTEST*/

