#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/moments.h"
#include "../../core/correlation.h"
#include "../../minimal/common.h"
#include "../tests.h"

#define DISCRETE_SWAP_X() \
      xdata = VECTOR_ELT(xx, i); \
      xptr = INTEGER(xdata); \
      llx = NLEVELS(xdata);

#define GAUSSIAN_SWAP_X() \
      xptr = REAL(VECTOR_ELT(xx, i)); \
      xm = c_mean(xptr, nobs); \
      xsd = c_sse(xptr, xm, nobs);

/* parametric tests for discrete variables. */
double ut_discrete(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df, test_e test) {

int i = 0, llx = 0, lly = NLEVELS(yy), *xptr = NULL, *yptr = INTEGER(yy);
double statistic = 0;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_chisqtest(xptr, llx, yptr, lly, nobs, df, test,
                    (test == MI) || (test == MI_ADF));
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shmi(xptr, llx, yptr, lly, nobs, TRUE);
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
double ut_gaustests_complete(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0;
double transform = 0, *xptr = NULL, *yptr = REAL(yy);
double xm = 0, ym = 0, xsd = 0, ysd = 0, statistic = 0;

  /* compute the degrees of freedom for correlation and mutual information. */
  *df = gaussian_df(test, nobs);

  if (*df < 1) {

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

      statistic *= 1 - cor_lambda(xptr, yptr, nobs, nobs, xm, ym, xsd, ysd,
                         statistic);
      statistic = 2 * nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = cor_zf_trans(statistic, *df);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_GAUSTESTS_COMPLETE*/

double ut_gaustests_with_missing(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0, ncomplete = 0;
double transform = 0, *xptr = NULL, *yptr = REAL(yy);
double statistic = 0, xm = 0, ym = 0, xsd = 0, ysd = 0;

  for (i = 0; i < ntests; i++) {

    /* compute the correlation. */
    xptr = REAL(VECTOR_ELT(xx, i));
    statistic = c_cor_with_missing(xptr, yptr, nobs, &xm, &ym, &xsd, &ysd,
                  &ncomplete);

    /* compute the degrees of freedom for correlation and mutual information. */
    *df = gaussian_df(test, ncomplete);

    if ((ncomplete == 0) || (*df < 1)) {

      /* if there are not enough degrees of freedom, including when there are
       * no complete observations, return perfect independence. */
      warning("trying to do an independence test with zero degrees of freedom.");

      *df = 0;
      statistic = 0;
      pvalue[i] = 1;

      continue;

    }/*THEN*/

    if (test == COR) {

      transform = cor_t_trans(statistic, *df);
      pvalue[i] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      statistic *= 1 - cor_lambda(xptr, yptr, nobs, ncomplete, xm, ym, xsd, ysd, statistic);
      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = cor_zf_trans(statistic, *df);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*THEN*/

  return statistic;

}/*UT_GAUSTESTS_WITH_MISSING*/

/* conditional linear Gaussian mutual information test. */
double ut_micg_complete(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
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
      statistic = c_chisqtest(xptr, llx, yptr, lly, nobs, df, MI, TRUE);
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
        statistic = 2 * nobs * c_micg(yptr, ym, ysd, xptr, llx, nobs, df);
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*THEN*/
      else {

        xptr = REAL(xdata);
        xm = c_mean(xptr, nobs);
        xsd = sqrt(c_sse(xptr, xm, nobs) / (nobs - 1));
        statistic = 2 * nobs * c_micg(xptr, xm, xsd, yptr, lly, nobs, df);
        pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

      }/*ELSE*/

    }/*THEN*/

  }/*FOR*/

  return statistic;

}/*UT_MICG_COMPLETE*/

double ut_micg_with_missing(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df) {

int i = 0, nc = 0, xtype = 0, ytype = TYPEOF(yy), lly = 0, disc_lvls = 0;
int *disc_ptr = NULL;
double statistic = 0, *cont_ptr = NULL;
void *yptr = NULL;
SEXP xdata;

  yptr = DATAPTR(yy);
  if (ytype == INTSXP)
    lly = NLEVELS(yy);

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xtype = TYPEOF(xdata);

    if ((ytype == INTSXP) && (xtype == INTSXP)) {

      /* if both nodes are discrete, the test reverts back to a discrete
       * mutual information test which handles missing data internally. */
      statistic = c_chisqtest(INTEGER(xdata), NLEVELS(xdata), yptr, lly,
                    nobs, df, MI, TRUE);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      /* if both nodes are continuous, the test reverts back to a Gaussian
       * mutual information test which handles missing data in estimating
       * the correlation */
      statistic = c_cor_with_missing(REAL(xdata), yptr, nobs, NULL, NULL,
                    NULL, NULL, &nc);
      *df = 1;
      statistic = 2 * nc * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else {

      if (xtype == INTSXP) {

        disc_ptr = INTEGER(xdata);
        disc_lvls = NLEVELS(xdata);
        cont_ptr = yptr;

      }/*THEN*/
      else {

        disc_ptr = yptr;
        disc_lvls = lly;
        cont_ptr = REAL(xdata);

      }/*ELSE*/

      statistic = c_micg_with_missing(cont_ptr, disc_ptr, disc_lvls,
                    nobs, df, &nc);
      statistic *= 2 * nc;
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*ELSE*/

  }/*FOR*/

  return statistic;

}/*UT_MICG_WITH MISSING*/

/* discrete permutation tests. */
double ut_dperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df, test_e type, int B, double a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy);
int llx = 0, lly = NLEVELS(yy);
double statistic = 0;
SEXP xdata;

  for (i = 0; i < ntests; i++) {

    DISCRETE_SWAP_X();
    c_mcarlo(xptr, llx, yptr, lly, nobs, B, &statistic, pvalue + i,
      a, type, df);

  }/*FOR*/

  return statistic;

}/*UT_DPERM*/

/* continuous permutation tests. */
double ut_gperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    test_e type, int B, double a, bool complete) {

int i = 0, k = 0, nc = nobs;
double *yptr = REAL(yy), *xptr = NULL, *yptr_complete = NULL, *xptr_complete = NULL;
double statistic = 0;

  if (!complete) {

    xptr_complete = Calloc1D(nobs, sizeof(double));
    yptr_complete = Calloc1D(nobs, sizeof(double));

  }/*THEN*/
  else {

    yptr_complete = yptr;

  }/*ELSE*/

  for (i = 0; i < ntests; i++) {

    xptr = REAL(VECTOR_ELT(xx, i));

    if (!complete) {

      for (k = 0, nc = 0; k < nobs; k++) {

        if (ISNAN(xptr[k]) || ISNAN(yptr[k]))
          continue;

        xptr_complete[nc] = xptr[k];
        yptr_complete[nc] = yptr[k];
        nc++;

      }/*FOR*/

    }/*THEN*/
    else {

      xptr_complete = xptr;

    }/*ELSE*/

    c_gauss_mcarlo(xptr_complete, yptr_complete, nc, B, pvalue + i,
      a, type, &statistic);

  }/*FOR*/

  if (!complete) {

    Free1D(xptr_complete);
    Free1D(yptr_complete);

  }/*THEN*/

  return statistic;

}/*UT_GPERM*/

