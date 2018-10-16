#include "include/rcore.h"
#include "include/sets.h"
#include "include/data.frame.h"
#include "include/tests.h"
#include "include/covariance.h"
#include "include/globals.h"
#include "include/blas.h"
#include "include/matrix.h"
#include "include/data.table.h"

/* parametric tests for discrete variables. */
static double ct_discrete(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e test) {

int i = 0, llx = 0, lly = NLEVELS(yy), llz = 0;
int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
double statistic = 0;
SEXP xdata, config;

  PROTECT(config = c_configurations(zz, TRUE, TRUE));
  zptr = INTEGER(config);
  llz = NLEVELS(config);

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xptr = INTEGER(xdata);
    llx = NLEVELS(xdata);

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, nobs, df, test,
                    (test == MI) || (test == MI_ADF));
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, nobs, df, TRUE);
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
static double ct_gaustests_complete(SEXP xx, SEXP yy, SEXP zz, double *pvalue,
    double *df, test_e test) {

int i = 0, ntests = length(xx);
double transform = 0, statistic = 0, lambda = 0;
gdata dt = { 0 };
covariance cov = { 0 }, basecov = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 2);
  dt.col[1] = REAL(yy);

  /* compute the degrees of freedom for correlation and mutual information. */
  if (test == COR)
    *df = dt.m.nobs - dt.m.ncols;
  else if ((test == MI_G) || (test == MI_G_SH))
    *df = 1;

  if (((test == COR) && (*df < 1)) || ((test == ZF) && (dt.m.nobs - dt.m.ncols < 2))) {

    /* if there are not enough degrees of freedom, return independence. */
    warning("trying to do a conditional independence test with zero degrees of freedom.");

    *df = 0;
    statistic = 0;
    for (i = 0; i < ntests; i++)
      pvalue[i] = 1;

    return statistic;

  }/*THEN*/

  /* allocate the covariance matrix. */
  cov = new_covariance(dt.m.ncols, TRUE);
  basecov = new_covariance(dt.m.ncols, TRUE);

  /* compute the mean values and the covariance matrix. */
  gdata_cache_means(&dt, 1);
  c_covmat(dt.col, dt.mean, dt.m.nobs, dt.m.ncols, cov, 1);
  if (ntests > 1)
    copy_covariance(&cov, &basecov);

  /* compute the partial correlation and the test statistic. */
  for (i = 0; i < ntests; i++) {

    /* extract and plug-in the i-th variable. */
    dt.col[0] = REAL(VECTOR_ELT(xx, i));
    /* update the corresponding mean in the cache. */
    dt.mean[0] = c_mean(dt.col[0], dt.m.nobs);
    /* update the covariance matrix. */
    if (ntests > 1)
      copy_covariance(&basecov, &cov);
    c_update_covmat(dt.col, dt.mean, 0, dt.m.nobs, dt.m.ncols, cov.mat);

    if (test == COR) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      transform = cor_t_trans(statistic, *df);
      pvalue[i] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = 2 * dt.m.nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(dt.col, dt.mean, cov, dt.m.nobs, NULL,
                 dt.m.nobs);
      covmat_shrink(cov, lambda);
      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = 2 * dt.m.nobs * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = cor_zf_trans(statistic, (double)(dt.m.nobs - dt.m.ncols));
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  FreeCOV(basecov);
  FreeCOV(cov);
  FreeGDT(dt, FALSE);

  return statistic;

}/*CT_GAUSTESTS_COMPLETE*/

static double ct_gaustests_with_missing(SEXP xx, SEXP yy, SEXP zz,
    double *pvalue, double *df, test_e test) {

int i = 0, j = 0, ncomplete = 0, ntests = length(xx);
double transform = 0, *mean = NULL, statistic = 0, lambda = 0;
short int *missing_yz = NULL, *missing_all = NULL;
gdata dt = { 0 };
covariance cov = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 2);
  dt.col[1] = REAL(yy);
  /* allocate the mean vector. */
  mean = Calloc1D(dt.m.ncols, sizeof(double));
  /* allocate the covariance matrix. */
  cov = new_covariance(dt.m.ncols, TRUE);

  /* allocate the missing values indicators. */
  missing_yz = Calloc1D(dt.m.nobs, sizeof(short int));
  if (test == MI_G_SH)
    missing_all = Calloc1D(dt.m.nobs, sizeof(short int));

  for (i = 0; i < dt.m.nobs; i++)
    for (j = 1; j < dt.m.ncols; j++)
      if (ISNAN(dt.col[j][i])) {

        missing_yz[i] = TRUE;
        break;

      }/*THEN*/

  /* compute the partial correlation and the test statistic. */
  for (i = 0; i < ntests; i++) {

    /* extract and plug in the i-th variable. */
    dt.col[0] = REAL(VECTOR_ELT(xx, i));
    /* compute the covariance matrix. */
    c_covmat_with_missing(dt.col, dt.m.nobs, dt.m.ncols, missing_yz,
      missing_all, mean, cov.mat, &ncomplete);

    /* compute the degrees of freedom for correlation and mutual information. */
    if (test == COR)
      *df = ncomplete - dt.m.ncols;
    else if ((test == MI_G) || (test == MI_G_SH))
      *df = 1;

    if ((ncomplete == 0) || ((test == COR) && (*df < 1)) ||
        ((test == ZF) && (ncomplete - dt.m.ncols < 2))) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      *df = 0;
      statistic = 0;
      pvalue[i] = 1;

      continue;

    }/*THEN*/

    if (test == COR) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      transform = cor_t_trans(statistic, *df);
      pvalue[i] = 2 * pt(fabs(transform), *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(dt.col, mean, cov, dt.m.nobs, missing_all,
                 ncomplete);
      covmat_shrink(cov, lambda);
      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
      statistic = cor_zf_trans(statistic, (double)(ncomplete - dt.m.ncols));
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  Free1D(mean);
  Free1D(missing_yz);
  if (test == MI_G_SH)
    Free1D(missing_all);
  FreeCOV(cov);
  FreeGDT(dt, FALSE);

  return statistic;

}/*CT_GAUSTESTS_WITH_MISSING*/

/* conditional linear Gaussian mutual information test. */
static double ct_micg_complete(SEXP xx, SEXP yy, SEXP zz, double *pvalue,
    double *df) {

int i = 0, *zptr = 0, ntests = length(xx);
int xtype = 0, ytype = TYPEOF(yy), llx = 0, lly = 0, llz = 0;
void *xptr = NULL, *yptr = NULL;
double statistic = 0;
SEXP xdata;
cgdata dt = { 0 };

  if (ytype == INTSXP) {

    /* cache the number of levels. */
    lly = NLEVELS(yy);
    yptr = INTEGER(yy);

  }/*THEN*/
  else {

    yptr = REAL(yy);

  }/*ELSE*/

  /* allocate and initialize a data table for the variables. */
  dt = cgdata_from_SEXP(zz, 1, 1);

  /* allocate vector for the configurations of the discrete parents. */
  if (dt.ndcols - 1 > 0) {

    zptr = Calloc1D(dt.m.nobs, sizeof(int));
    c_fast_config(dt.dcol + 1, dt.m.nobs, dt.ndcols - 1, dt.nlvl + 1,
      zptr, &llz, 1);

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

      if (dt.ngcols - 1 > 0) {

        /* need to reverse conditioning to actually compute the test. */
        statistic = 2 * dt.m.nobs * dt.m.nobs *
                      c_cmicg_unroll(xptr, llx, yptr, lly, zptr, llz,
                               dt.gcol + 1, dt.ngcols - 1, df, dt.m.nobs);

      }/*THEN*/
      else {

        /* the test reverts back to a discrete mutual information test. */
        statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, dt.m.nobs, df,
                      MI, TRUE);

      }/*ELSE*/

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      dt.gcol[0] = xptr;
      statistic = 2 * dt.m.nobs * c_cmicg(yptr, dt.gcol, dt.ngcols, NULL,
                        0, zptr, llz, dt.nlvl, dt.m.nobs, df);

    }/*THEN*/
    else if ((ytype == INTSXP) && (xtype == REALSXP)) {

      dt.dcol[0] = yptr;
      dt.nlvl[0] = lly;
      statistic = 2 * dt.m.nobs * c_cmicg(xptr, dt.gcol + 1, dt.ngcols - 1,
                        dt.dcol, dt.ndcols, zptr, llz, dt.nlvl,
                        dt.m.nobs, df);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == INTSXP)) {

      dt.dcol[0] = xptr;
      dt.nlvl[0] = llx;
      statistic = 2 * dt.m.nobs * c_cmicg(yptr, dt.gcol + 1, dt.ngcols - 1,
                        dt.dcol, dt.ndcols, zptr, llz, dt.nlvl,
                        dt.m.nobs, df);

    }/*ELSE*/

    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  Free1D(zptr);
  FreeCGDT(dt, FALSE);

  return statistic;

}/*CT_MICG_COMPLETE*/

static double ct_micg_with_missing(SEXP xx, SEXP yy, SEXP zz, double *pvalue,
    double *df) {

int xtype = 0, ytype = TYPEOF(yy), llx = 0, lly = 0, llz = 0, ntests = length(xx);
int i = 0, j = 0, k = 0, l = 0, nc = 0, *zptr = 0, **dp_complete = NULL;
void *xptr = NULL, *xptr_complete = NULL, *yptr = NULL, *yptr_complete = NULL;
double **gp_complete = NULL;
double statistic = 0;
short int *missing_yz = NULL, *missing_all = NULL;
SEXP xdata;
cgdata dt = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = cgdata_from_SEXP(zz, 1, 1);

  if (ytype == INTSXP) {

    /* cache the number of levels. */
    lly = NLEVELS(yy);
    yptr = INTEGER(yy);
    yptr_complete = Calloc1D(dt.m.nobs, sizeof(int));

  }/*THEN*/
  else {

    yptr = REAL(yy);
    yptr_complete = Calloc1D(dt.m.nobs, sizeof(double));

  }/*ELSE*/

  /* extract the conditioning variables and cache their types. */
  dp_complete = (int **)Calloc2D(dt.ndcols, dt.m.nobs, sizeof(int));
  gp_complete = (double **)Calloc2D(dt.ngcols, dt.m.nobs, sizeof(double));

  /* compute the missing values indicators. */
  missing_yz = Calloc1D(dt.m.nobs, sizeof(short int));
  missing_all = Calloc1D(dt.m.nobs, sizeof(short int));

  for (i = 0; i < dt.m.nobs; i++) {

    if (ytype == INTSXP) {

      if (((int *)yptr)[i] == NA_INTEGER) {

        missing_yz[i] = TRUE;
        continue;

      }/*THEN*/

    }/*THEN*/
    else {

      if (ISNAN(((double *)yptr)[i])) {

        missing_yz[i] = TRUE;
        continue;

      }/*THEN*/

    }/*ELSE*/

    for (j = 0; j < dt.ndcols - 1; j++)
      if (dt.dcol[1 + j][i] == NA_INTEGER) {

        missing_yz[i] = TRUE;
        continue;

      }/*THEN*/

    for (j = 0; j < dt.ngcols - 1; j++)
      if (ISNAN(dt.gcol[1 + j][i])) {

        missing_yz[i] = TRUE;
        continue;

      }/*THEN*/

  }/*FOR*/

  /* allocate vector for the configurations of the discrete parents. */
  if (dt.ndcols - 1 > 0)
    zptr = Calloc1D(dt.m.nobs, sizeof(int));

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xtype = TYPEOF(xdata);

    if (xtype == INTSXP) {

      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      for (j = 0; j < dt.m.nobs; j++)
        missing_all[j] = missing_yz[j] || (((int *)xptr)[j] == NA_INTEGER);

      xptr_complete = Calloc1D(dt.m.nobs, sizeof(int));
      for (j = 0, k = 0; j < dt.m.nobs; j++)
        if (!missing_all[j])
          ((int *)xptr_complete)[k++] = ((int *)xptr)[j];

    }/*THEN*/
    else {

      xptr = REAL(xdata);
      for (j = 0; j < dt.m.nobs; j++)
        missing_all[j] = missing_yz[j] || ISNAN(((double *)xptr)[j]);

      xptr_complete = Calloc1D(dt.m.nobs, sizeof(double));
      for (j = 0, k = 0; j < dt.m.nobs; j++)
        if (!missing_all[j])
          ((double *)xptr_complete)[k++] = ((double *)xptr)[j];

    }/*ELSE*/

    /* subset the complete data. */
    for (j = 0; j < dt.m.nobs; j++)
      nc += !missing_all[j];

    for (j = 0; j < dt.ngcols - 1; j++)
      for (l = 0, k = 0; l < dt.m.nobs; l++)
        if (!missing_all[l])
          gp_complete[1 + j][k++] = dt.gcol[1 + j][l];

    if (dt.ndcols - 1 > 0) {

      for (j = 0; j < dt.ndcols - 1; j++)
        for (l = 0, k = 0; l < dt.m.nobs; l++)
          if (!missing_all[l])
            dp_complete[1 + j][k++] = dt.dcol[1 + j][l];

      c_fast_config(dp_complete + 1, nc, dt.ndcols - 1, dt.nlvl + 1, zptr,
        &llz, 1);

    }/*THEN*/

    if (ytype == INTSXP) {

      for (j = 0, k = 0; j < dt.m.nobs; j++)
        if (!missing_all[j])
          ((int *)yptr_complete)[k++] = ((int *)yptr)[j];

    }/*THEN*/
    else {

      for (j = 0, k = 0; j < dt.m.nobs; j++)
        if (!missing_all[j])
          ((double *)yptr_complete)[k++] = ((double *)yptr)[j];

    }/*ELSE*/

    if (nc == 0) {

      statistic = 0;
      pvalue[i] = 1;
      Free1D(xptr_complete);
      continue;

    }/*THEN*/

    if ((ytype == INTSXP) && (xtype == INTSXP)) {

      if (dt.ngcols - 1 > 0) {

        /* need to reverse conditioning to actually compute the test. */
        statistic = 2 * nc * nc *
                      c_cmicg_unroll(xptr_complete, llx, yptr_complete, lly,
                        zptr, llz, gp_complete + 1, dt.ngcols - 1, df, nc);

      }/*THEN*/
      else {

        /* the test reverts back to a discrete mutual information test. */
        statistic = c_cchisqtest(xptr_complete, llx, yptr_complete, lly, zptr,
                      llz, nc, df, MI, TRUE);

      }/*ELSE*/

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == REALSXP)) {

      memcpy(gp_complete[0], xptr_complete, nc * sizeof(double));
      statistic = 2 * nc * c_cmicg(yptr_complete, gp_complete, dt.ngcols, NULL,
                               0, zptr, llz, dt.nlvl, nc, df);

    }/*THEN*/
    else if ((ytype == INTSXP) && (xtype == REALSXP)) {

      memcpy(dp_complete[0], yptr_complete, nc * sizeof(int));
      dt.nlvl[0] = lly;
      statistic = 2 * nc * c_cmicg(xptr_complete, gp_complete + 1,
                             dt.ngcols - 1, dp_complete, dt.ndcols, zptr, llz,
                             dt.nlvl, nc, df);

    }/*THEN*/
    else if ((ytype == REALSXP) && (xtype == INTSXP)) {

      memcpy(dp_complete[0], xptr_complete, nc * sizeof(int));
      dt.nlvl[0] = llx;
      statistic = 2 * nc * c_cmicg(yptr_complete, gp_complete + 1,
                             dt.ngcols - 1, dp_complete, dt.ndcols, zptr,
                             llz, dt.nlvl, nc, df);

    }/*ELSE*/

    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    Free1D(xptr_complete);

  }/*FOR*/

  Free1D(zptr);
  Free2D(dp_complete, dt.ndcols);
  Free2D(gp_complete, dt.ngcols);
  Free1D(yptr_complete);
  Free1D(missing_yz);
  Free1D(missing_all);
  FreeCGDT(dt, FALSE);

  return statistic;

}/*CT_MICG_WITH_MISSING*/

/* discrete permutation tests. */
static double ct_dperm(SEXP xx, SEXP yy, SEXP zz, int nobs, int ntests,
    double *pvalue, double *df, test_e type, int B, double a) {

int i = 0, *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
int llx = 0, lly = NLEVELS(yy), llz = 0;
double statistic = 0;
SEXP xdata, config;

  PROTECT(config = c_configurations(zz, TRUE, TRUE));
  zptr = INTEGER(config);
  llz = NLEVELS(config);

  for (i = 0; i < ntests; i++) {

    xdata = VECTOR_ELT(xx, i);
    xptr = INTEGER(xdata);
    llx = NLEVELS(xdata);
    statistic = 0;
    c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, B, &statistic,
      pvalue + i, a, type, df);

  }/*FOR*/

  UNPROTECT(1);

  return statistic;

}/*CT_DPERM*/

/* continuous permutation tests. */
static double ct_gperm(SEXP xx, SEXP yy, SEXP zz, double *pvalue, double *df,
    test_e type, int B, double a, int complete) {

int i = 0, j = 0, k = 0, nc = 0, ntests = length(xx);
double **complete_column = NULL, *yptr = REAL(yy), statistic = 0;
short int *missing_yz = NULL;
gdata dt = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 2);
  dt.col[1] = REAL(yy);

  if (!complete) {

    /* allocate the missing values indicators. */
    missing_yz = Calloc1D(dt.m.nobs, sizeof(short int));

    for (i = 0; i < dt.m.nobs; i++)
      for (j = 1; j < dt.m.ncols; j++)
        if (ISNAN(dt.col[j][i])) {

          missing_yz[i] = TRUE;
          break;

        }/*THEN*/

    complete_column = (double **) Calloc2D(dt.m.ncols, dt.m.nobs, sizeof(double));

  }/*THEN*/
  else {

    complete_column = dt.col;

  }/*ELSE*/

  for (i = 0; i < ntests; i++) {

    /* swap the first column and restore the second, which is that undergoing
     * permutation (backward compatibility from set random seed). */
    dt.col[0] = REAL(VECTOR_ELT(xx, i));
    dt.col[1] = yptr;

    if (!complete) {

      for (k = 0, nc = 0; k < dt.m.nobs; k++) {

        if (missing_yz[k] || ISNAN(dt.col[0][k]))
          continue;

        for (j = 0; j < dt.m.ncols; j++)
          complete_column[j][nc] = dt.col[j][k];
        nc++;

      }/*FOR*/

    }/*THEN*/
    else {

      nc = dt.m.nobs;

    }/*ELSE*/

    statistic = 0;
    c_gauss_cmcarlo(complete_column, dt.m.ncols, nc, 0, 1, B, &statistic,
      pvalue + i, a, type);

  }/*FOR*/

  if (!complete) {

    Free1D(missing_yz);
    Free2D(complete_column, dt.m.ncols);

  }/*THEN*/

  FreeGDT(dt, FALSE);

  return statistic;

}/*CT_GPERM*/

/* conditional independence tests. */
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning, SEXP complete) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
SEXP xx, yy, zz, cc, result;

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

  /* extract the missing values indicators. */
  PROTECT(cc = subset_by_name(complete, 3, y, x, sx));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    statistic = ct_discrete(xx, yy, zz, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    if (all_equal(cc, TRUESEXP))
      statistic = ct_gaustests_complete(xx, yy, zz, pvalue, &df, test_type);
    else
      statistic = ct_gaustests_with_missing(xx, yy, zz, pvalue, &df, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    if (all_equal(cc, TRUESEXP))
      statistic = ct_micg_complete(xx, yy, zz, pvalue, &df);
    else
      statistic = ct_micg_with_missing(xx, yy, zz, pvalue, &df);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    statistic = ct_dperm(xx, yy, zz, nobs, ntests, pvalue, &df, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    statistic = ct_gperm(xx, yy, zz, pvalue, &df, test_type, INT(B),
                  IS_SMC(test_type) ? NUM(alpha) : 1, all_equal(cc, TRUESEXP));

  }/*THEN*/

  UNPROTECT(5);

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

