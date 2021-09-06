#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/data.table.h"
#include "include/matrix.h"

static void rrd_disc_message(meta m, const char *x, int offset, const char *y,
    double p, double alpha) {

  Rprintf("    > node %s is %s %s given ", x,
    (p <= alpha ? "dependent on" : "independent from"), y);
  for (int j = offset; j < m.ncols; j++)
    Rprintf("%s ", m.names[j]);
  Rprintf("(p-value: %g).\n", p);

}/*RRD_DISC_MESSAGE*/

static void rrd_gauss_message(gdata dt, int i, double p, double alpha) {

  Rprintf("    > node %s is %s %s given ", dt.m.names[0],
    (p <= alpha ? "dependent on" : "independent from"), dt.m.names[i]);
  for (int j = 1; j < dt.m.ncols; j++)
    if (j != i)
      Rprintf("%s ", dt.m.names[j]);
  Rprintf("(p-value: %g).\n", p);

}/*RRD_GAUSS_MESSAGE*/

/* parametric tests for discrete variables. */
static void rrd_discrete(ddata dtx, ddata dtz, test_e test, double *pvalue,
    double alpha, bool debugging) {

int i = 0, cur = 0, valid = dtz.m.ncols, lly = 0, llz = 0;
int *yptr = NULL, *zptr = NULL;
double statistic = 0, df = 0;
ddata sub = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dtz.m.nobs, dtz.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dtz.m.nobs, sizeof(int));

  for (i = 0; i < dtz.m.ncols; i++) {

    /* nothing to do if there is just one variable left. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dtz.m.flag[i].fixed)
      continue;

    /* extract the variable to test. */
    yptr = dtz.col[i];
    lly = dtz.nlvl[i];
    /* drop the variable to test from the conditioning set. */
    dtz.m.flag[i].drop = TRUE;
    ddata_drop_flagged(&dtz, &sub);
    /* construct the configurations of the conditioning variables. */
    c_fast_config(sub.col, sub.m.nobs, sub.m.ncols, sub.nlvl, zptr, &llz, 1);

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(dtx.col[0], dtx.nlvl[0], yptr, lly, zptr, llz,
                    sub.m.nobs, &df, test, (test == MI) || (test == MI_ADF));
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shcmi(dtx.col[0], dtx.nlvl[0], yptr, lly, zptr, llz,
                    sub.m.nobs, &df, TRUE);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_cjt(dtx.col[0], dtx.nlvl[0], yptr, lly, zptr, llz,
                    sub.m.nobs);
      pvalue[cur] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    if (debugging)
      rrd_disc_message(sub.m, dtx.m.names[0], 0, dtz.m.names[i], pvalue[cur],
        alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dtz.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(sub);

}/*RRD_DISCRETE*/

/* parametric tests for gaussian variables (and complete data). */
static void rrd_gaustests_complete(gdata dt, test_e test, double *pvalue,
    double alpha, bool debugging) {

int i = 1, cur = 0, df = 0, valid = dt.m.ncols - 1, t = 0, run_svd = TRUE;
double statistic = 0, lambda = 0;
gdata sub = { 0 };
covariance cov = { 0 }, backup = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.m.names = Calloc1D(dt.m.ncols, sizeof(char *));
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));

  /* allocate the covariance matrix. */
  cov = new_covariance(dt.m.ncols, TRUE);
  backup = new_covariance(dt.m.ncols, TRUE);
  c_covmat(dt.col, dt.mean, dt.m.nobs, dt.m.ncols, cov, 0);

  /* the counter starts at 1 because the first column is the target variable,
   * swapping it does not make sense and would just add a duplicate, redundant
   * test. */
  for (i = 1; i < dt.m.ncols; i++) {

    /* nothing to do if there is just one variable. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* drop the variables that have already been found to be independent. */
    gdata_drop_flagged(&dt, &sub);
    t = i - (dt.m.ncols - sub.m.ncols);

    /* compute the degrees of freedom for correlation and mutual information. */
    df = gaussian_cdf(test, sub.m.nobs, sub.m.ncols - 2);

    if (df < 1) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      pvalue[cur++] = 1;

      if (debugging)
        rrd_gauss_message(sub, t, 1, alpha);

      continue;

    }/*THEN*/

    /* backup a copy of the covariance matrix before messing with it. */
    copy_covariance(&cov, &backup);

    if (test == COR) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = cor_t_trans(statistic, (double)df);
      pvalue[cur] = 2 * pt(fabs(statistic), df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(sub.col, sub.mean, cov, sub.m.nobs,
                 NULL, sub.m.nobs);
      covmat_shrink(cov, lambda);
      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = cor_zf_trans(statistic, df);
      pvalue[cur] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    if (debugging)
      rrd_gauss_message(sub, t, pvalue[cur], alpha);

    /* restore the covariance matrix before subsetting it. */
    copy_covariance(&backup, &cov);

    /* remove ifrom the covariance matrix and make sure to recompute SVD. */
    if (pvalue[cur++] > alpha) {

      valid--;
      dt.m.flag[i].drop = TRUE;
      covariance_drop_variable(&cov, &cov, t);
      run_svd = TRUE;

    }/*THEN*/
    else {

      run_svd = FALSE;

    }/*ELSE*/

  }/*WHILE*/

  FreeGDT(sub);
  FreeCOV(backup);
  FreeCOV(cov);

}/*RRD_GAUSTESTS_COMPLETE*/

/* parametric tests for gaussian variables (and incomplete data). */
static void rrd_gaustests_with_missing(gdata dt, test_e test, double *pvalue,
    double alpha, bool debugging) {

int i = 1, cur = 0, df = 0, valid = dt.m.ncols - 1, t = 0, run_svd = TRUE;
int ncomplete = 0;
double statistic = 0, lambda = 0, *mean = NULL;
bool *missing = NULL;
gdata sub = { 0 };
covariance cov = { 0 }, backup = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.m.names = Calloc1D(dt.m.ncols, sizeof(char *));

  /* allocate the covariance matrix. */
  cov = new_covariance(dt.m.ncols, TRUE);
  backup = new_covariance(dt.m.ncols, TRUE);

  missing = Calloc1D(dt.m.nobs, sizeof(bool));
  mean = Calloc1D(dt.m.ncols, sizeof(double));
  c_covmat_with_missing(dt.col, dt.m.nobs, dt.m.ncols, NULL, missing, mean,
    cov.mat, &ncomplete);

  /* the counter starts at 1 because the first column is the target variable,
   * swapping it does not make sense and would just add a duplicate, redundant
   * test. */
  for (i = 1; i < dt.m.ncols; i++) {

    /* nothing to do if there is just one variable. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* drop the variables that have already been found to be independent. */
    gdata_drop_flagged(&dt, &sub);
    t = i - (dt.m.ncols - sub.m.ncols);

    /* compute the degrees of freedom for correlation and mutual information. */
    df = gaussian_cdf(test, ncomplete, sub.m.ncols - 2);

    if ((ncomplete == 0) || (df < 1)) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      pvalue[cur++] = 1;

      if (debugging)
        rrd_gauss_message(sub, t, 1, alpha);

      continue;

    }/*THEN*/

    /* backup a copy of the covariance matrix before messing with it. */
    copy_covariance(&cov, &backup);

    if (test == COR) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = cor_t_trans(statistic, (double)df);
      pvalue[cur] = 2 * pt(fabs(statistic), df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_G_SH) {

      lambda = covmat_lambda(sub.col, mean, cov, sub.m.nobs, missing, ncomplete);
      covmat_shrink(cov, lambda);
      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = 2 * ncomplete * cor_mi_trans(statistic);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == ZF) {

      statistic = c_fast_pcor(cov, 0, t, NULL, run_svd);
      statistic = cor_zf_trans(statistic, df);
      pvalue[cur] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    if (debugging)
      rrd_gauss_message(sub, t, pvalue[cur], alpha);

    /* restore the covariance matrix before subsetting it. */
    copy_covariance(&backup, &cov);

    /* remove ifrom the covariance matrix and make sure to recompute SVD. */
    if (pvalue[cur++] > alpha) {

      valid--;
      dt.m.flag[i].drop = TRUE;
      gdata_drop_flagged(&dt, &sub);
      c_covmat_with_missing(sub.col, sub.m.nobs, sub.m.ncols, NULL, missing,
        mean, cov.mat, &ncomplete);

      run_svd = TRUE;

    }/*THEN*/
    else {

      run_svd = FALSE;

    }/*ELSE*/

  }/*WHILE*/

  FreeGDT(sub);
  FreeCOV(backup);
  FreeCOV(cov);
  Free1D(missing);
  Free1D(mean);

}/*RRD_GAUSTESTS_WITH_MISSING*/

/* conditional linear Gaussian test (for complete data). */
double rrd_micg_chisq(cgdata dtx, cgdata dty, cgdata sub, int *zptr, int llz,
    double *df, bool copy) {

double statistic = 0;

  if (dtx.m.flag[0].discrete && dty.m.flag[0].discrete) {

    /* check whether the conditioning set is valid. */
    if (sub.ngcols > 1) {

      /* need to reverse conditioning to actually compute the test. */
      statistic = 2 * sub.m.nobs * sub.m.nobs *
                    c_cmicg_unroll(dtx.dcol[0], dtx.nlvl[0], dty.dcol[0],
                      dty.nlvl[0], zptr, llz, sub.gcol + 1, sub.ngcols - 1, df,
                      sub.m.nobs);

    }/*THEN*/
    else {

      /* if both nodes are discrete, the test reverts back to a discrete
       * mutual information test. */
      statistic = c_cchisqtest(dtx.dcol[0], dtx.nlvl[0], dty.dcol[0],
                    dty.nlvl[0], zptr, llz, sub.m.nobs, df, MI, TRUE);

    }/*ELSE*/

  }/*THEN*/
  else if (dtx.m.flag[0].gaussian && dty.m.flag[0].gaussian) {

    if (copy)
      memcpy(sub.gcol[0], dtx.gcol[0], sub.m.nobs * sizeof(double));
    else
      sub.gcol[0] = dtx.gcol[0];
    statistic = 2 * sub.m.nobs *
                  c_cmicg(dty.gcol[0], sub.gcol, sub.ngcols, NULL, 0, zptr, llz,
                    0, sub.m.nobs, df);

  }/*THEN*/
  else if (dtx.m.flag[0].gaussian && dty.m.flag[0].discrete) {

    if (copy)
      memcpy(sub.dcol[0], dty.dcol[0], sub.m.nobs * sizeof(int));
    else
      sub.dcol[0] = dty.dcol[0];
    sub.nlvl[0] = dty.nlvl[0];
    statistic = 2 * sub.m.nobs *
                  c_cmicg(dtx.gcol[0], sub.gcol + 1, sub.ngcols - 1, sub.dcol,
                           sub.ndcols, zptr, llz, sub.nlvl, sub.m.nobs, df);

  }/*THEN*/
  else if (dtx.m.flag[0].discrete && dty.m.flag[0].gaussian) {

    if (copy)
      memcpy(sub.dcol[0], dtx.dcol[0], sub.m.nobs * sizeof(int));
    else
      sub.dcol[0] = dtx.dcol[0];
    sub.nlvl[0] = dtx.nlvl[0];
    statistic = 2 * sub.m.nobs *
                  c_cmicg(dty.gcol[0], sub.gcol + 1, sub.ngcols - 1, sub.dcol,
                    sub.ndcols, zptr, llz, sub.nlvl, sub.m.nobs, df);

  }/*THEN*/

  return statistic;

}/*RRD_MICG_CHISQ*/

static void rrd_micg_complete(cgdata dtx, cgdata dtz, test_e test,
    double *pvalue, double alpha, bool debugging) {

int i = 0;
int cur = 0, valid = dtz.m.ncols - 2, llz = 0;
int *configurations = NULL, *zptr = NULL;
double statistic = 0, df = 0;
cgdata sub = { 0 }, dty = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);
  configurations = Calloc1D(dtz.m.nobs, sizeof(int));
  dty = empty_cgdata(dtz.m.nobs, 1, 1);

  for (i = 2; i < dtz.m.ncols; i++) {

    /* nothing to do if there are no more valid nodes.*/
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dtz.m.flag[i].fixed)
      continue;

    /* drop the variable to test from the conditioning set (along with variables
     * that have already been found to be independent). */
    dtz.m.flag[i].drop = TRUE;
    cgdata_drop_flagged(&dtz, &sub);
    /* extract the variable to test. */
    cgdata_subset_columns(&dtz, &dty, &i, 1);

    /* if there are discrete conditioning variables, compute their
     * configurations. */
    if (sub.ndcols  > 1) {

      zptr = configurations;
      c_fast_config(sub.dcol + 1, sub.m.nobs, sub.ndcols - 1, sub.nlvl + 1,
        zptr, &llz, 1);

    }/*THEN*/
    else {

      zptr = NULL;
      llz = 0;

    }/*ELSE*/

    statistic = rrd_micg_chisq(dtx, dty, sub, zptr, llz, &df, FALSE);

    pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    if (debugging)
      rrd_disc_message(sub.m, dtx.m.names[0], 2, dtz.m.names[i], pvalue[cur],
        alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dtz.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(configurations);
  FreeCGDT(sub);
  FreeCGDT(dty);

}/*RRD_MICG_COMPLETE*/

/* conditional linear Gaussian test (for incomplete data). */
static void rrd_micg_with_missing(cgdata dtx, cgdata dtz, test_e test,
    double *pvalue, double alpha, bool debugging) {

int i = 0;
int cur = 0, valid = dtz.m.ncols - 2, llz = 0;
int *configurations = NULL, *zptr = NULL;
double statistic = 0, df = 0;
bool *missing_x = NULL, *missing_all = NULL;
cgdata dtx_complete = { 0 }, dty = { 0 }, dty_complete = { 0 };
cgdata sub = { 0 }, sub_complete = { 0 };
meta sub_meta = { 0 }, dty_meta = { 0 }, dtx_complete_meta = { 0 };
meta dty_complete_meta = { 0 }, sub_complete_meta = { 0 };
int *sub_map = NULL, *dty_map = NULL, *dtx_complete_map = NULL;
int *dty_complete_map = NULL, *sub_complete_map = NULL;

  /* allocate a second data table to hold the conditioning variables, and back
   * their metadata up to restore them later. */
  dty = empty_cgdata(dtz.m.nobs, 1, 1);
  dty_meta.flag = Calloc1D(dty.m.ncols, sizeof(flags));
  meta_copy(&(dty.m), &dty_meta);
  dty_map = Calloc1D(dty.m.ncols, sizeof(int));
  memcpy(dty_map, dty.map, dty.m.ncols * sizeof(int));

  sub = empty_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);
  sub_meta.flag = Calloc1D(sub.m.ncols, sizeof(flags));
  meta_copy(&(sub.m), &sub_meta);
  sub_map = Calloc1D(sub.m.ncols, sizeof(int));
  memcpy(sub_map, sub.map, sub.m.ncols * sizeof(int));

  sub_complete = new_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);
  sub_complete_meta.flag = Calloc1D(sub_complete.m.ncols, sizeof(flags));
  meta_copy(&(sub_complete.m), &sub_complete_meta);
  sub_complete_map = Calloc1D(sub_complete.m.ncols, sizeof(int));
  memcpy(sub_complete_map, sub_complete.map, sub_complete.m.ncols * sizeof(int));

  dtx_complete = new_cgdata(dtx.m.nobs, dtx.ndcols, dtx.ngcols);
  dtx_complete_meta.flag = Calloc1D(dtx_complete.m.ncols, sizeof(flags));
  meta_copy(&(dtx_complete.m), &dtx_complete_meta);
  dtx_complete_map = Calloc1D(dtx_complete.m.ncols, sizeof(int));
  memcpy(dtx_complete_map, dtx_complete.map, dtx_complete.m.ncols * sizeof(int));

  dty_complete = new_cgdata(dty.m.nobs, dty.ndcols, dty.ngcols);
  dty_complete_meta.flag = Calloc1D(dty_complete.m.ncols, sizeof(flags));
  meta_copy(&(dty_complete.m), &dty_complete_meta);
  dty_complete_map = Calloc1D(dty_complete.m.ncols, sizeof(int));
  memcpy(dty_complete_map, dty_complete.map, dty_complete.m.ncols * sizeof(int));

  configurations = Calloc1D(dtz.m.nobs, sizeof(int));

  /* allocate the missingness indicators. */
  missing_x = Calloc1D(dtz.m.nobs, sizeof(bool));
  missing_all = Calloc1D(dtz.m.nobs, sizeof(bool));
  cgdata_incomplete_cases(&dtx, missing_x, 0, 0);

  for (i = 2; i < dtz.m.ncols; i++) {

    /* nothing to do if there are no more valid nodes.*/
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dtz.m.flag[i].fixed)
      continue;

    /* drop the variable to test from the conditioning set (along with variables
     * that have already been found to be independent). */
    dtz.m.flag[i].drop = TRUE;
    cgdata_drop_flagged(&dtz, &sub);

    /* extract the variable to test. */
    cgdata_subset_columns(&dtz, &dty, &i, 1);
    /* find out which observations are missing... */
    memcpy(missing_all, missing_x, dtz.m.nobs * sizeof(bool));
    cgdata_incomplete_cases(&sub, missing_all, 1, 1);
    cgdata_incomplete_cases(&dty, missing_all, 0, 0);
    /* ... and subset the complete data. */
    cgdata_subsample_by_logical(&sub, &sub_complete, missing_all, 1, 1);
    cgdata_subsample_by_logical(&dtx, &dtx_complete, missing_all, 0, 0);
    cgdata_subsample_by_logical(&dty, &dty_complete, missing_all, 0, 0);

    /* assume independence and return if there are no complete observations. */
    if (sub_complete.m.nobs == 0) {

      pvalue[cur] = 1;
      goto end;

    }/*THEN*/

    /* if there are discrete conditioning variables, compute their
     * configurations. */
    if (sub_complete.ndcols  > 1) {

      zptr = configurations;
      c_fast_config(sub_complete.dcol + 1, sub_complete.m.nobs,
        sub_complete.ndcols - 1, sub_complete.nlvl + 1, zptr, &llz, 1);

    }/*THEN*/
    else {

      zptr = NULL;
      llz = 0;

    }/*ELSE*/

    statistic = rrd_micg_chisq(dtx_complete, dty_complete, sub_complete, zptr,
                  llz, &df, TRUE);
    pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

end:

    if (debugging)
      rrd_disc_message(sub_complete.m, dtx_complete.m.names[0], 2,
        dtz.m.names[i], pvalue[cur], alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dtz.m.flag[i].drop = FALSE;

  }/*FOR*/

  /* restore data tables to their original state to ensure memory is freed.*/
  meta_copy(&dty_meta, &(dty.m));
  memcpy(dty.map, dty_map, dty.m.ncols * sizeof(int));
  meta_copy(&sub_meta, &(sub.m));
  memcpy(sub.map, sub_map, sub.m.ncols * sizeof(int));
  meta_copy(&sub_complete_meta, &(sub_complete.m));
  memcpy(sub_complete.map, sub_complete_map, sub_complete.m.ncols * sizeof(int));
  meta_copy(&dtx_complete_meta, &(dtx_complete.m));
  memcpy(dtx_complete.map, dtx_complete_map, dtx_complete.m.ncols * sizeof(int));
  meta_copy(&dty_complete_meta, &(dty_complete.m));
  memcpy(dty_complete.map, dty_complete_map, dty_complete.m.ncols * sizeof(int));

  Free1D(dty_map);
  Free1D(sub_map);
  Free1D(sub_complete_map);
  Free1D(dtx_complete_map);
  Free1D(dty_complete_map);
  FreeMETA(&dty_meta);
  FreeMETA(&sub_meta);
  FreeMETA(&sub_complete_meta);
  FreeMETA(&dtx_complete_meta);
  FreeMETA(&dty_complete_meta);
  Free1D(missing_x);
  Free1D(missing_all);
  Free1D(configurations);
  FreeCGDT(sub_complete);
  FreeCGDT(dtx_complete);
  FreeCGDT(dty_complete);
  FreeCGDT(dty);
  FreeCGDT(sub);

}/*RRD_MICG_WITH_MISSING*/

/* discrete permutation tests. */
static void rrd_dperm(ddata dtx, ddata dtz, test_e test, double *pvalue,
    double alpha, int nperms, double threshold, bool debugging) {

int i = 0, cur = 0, valid = dtz.m.ncols, lly = 0, llz = 0;
int *yptr = NULL, *zptr = NULL;
double statistic = 0, df = 0;
ddata sub = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dtz.m.nobs, dtz.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dtz.m.nobs, sizeof(int));

  for (i = 0; i < dtz.m.ncols; i++) {

    /* nothing to do if there is just one variable. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dtz.m.flag[i].fixed)
      continue;

    /* extract the variable to test. */
    yptr = dtz.col[i];
    lly = dtz.nlvl[i];
    /* drop the variable to test from the conditioning set. */
    dtz.m.flag[i].drop = TRUE;
    ddata_drop_flagged(&dtz, &sub);
    /* construct the configurations of the conditioning variables. */
    c_fast_config(sub.col, sub.m.nobs, sub.m.ncols, sub.nlvl, zptr, &llz, 1);

    c_cmcarlo(dtx.col[0], dtx.nlvl[0], yptr, lly, zptr, llz, sub.m.nobs, nperms,
      &statistic, pvalue + cur, threshold, test, &df);

    if (debugging)
      rrd_disc_message(sub.m, dtx.m.names[0], 0, dtz.m.names[i], pvalue[cur],
        alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dtz.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(sub);

}/*RRD_DPERM*/

/* continuous permutation tests. */
static void rrd_gperm(gdata dt, test_e test, double *pvalue, double alpha,
    int nperms, double threshold, bool complete, bool debugging) {

int i = 0, cur = 0, valid = dt.m.ncols - 1, t = 0;
double statistic = 0;
bool *missing = NULL;
gdata sub = { 0 }, sub_complete = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.m.names = Calloc1D(dt.m.ncols, sizeof(char *));

  if (!complete) {

    missing = Calloc1D(dt.m.nobs, sizeof(bool));
    sub_complete = new_gdata(dt.m.nobs, dt.m.ncols);
    sub_complete.m.names = Calloc1D(dt.m.ncols, sizeof(char *));

  }/*THEN*/

  /* the counter starts at 1 because the first column is the target variable,
   * swapping it does not make sense and would just add a duplicate, redundant
   * test. */
  for (i = 1; i < dt.m.ncols; i++) {

    /* nothing to do if there is just one variable. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* drop the variables that have already been found to be independent. */
    gdata_drop_flagged(&dt, &sub);
    t = i - (dt.m.ncols - sub.m.ncols);

    if (!complete) {

      gdata_incomplete_cases(&sub, missing, 0);
      gdata_subsample_by_logical(&sub, &sub_complete, missing, 0);
      c_gauss_cmcarlo(sub_complete.col, sub_complete.m.ncols,
        sub_complete.m.nobs, 0, t, nperms, &statistic, pvalue + cur, threshold,
        test);

    }/*THEN*/
    else {

      c_gauss_cmcarlo(sub.col, sub.m.ncols, sub.m.nobs, 0, t, nperms,
        &statistic, pvalue + cur, threshold, test);

    }/*ELSE*/

    if (debugging)
      rrd_gauss_message(sub, t, pvalue[cur], alpha);

    /* remove from the set of valid candidates for the conditioning set. */
    if (pvalue[cur++] > alpha) {

      valid--;
      dt.m.flag[i].drop = TRUE;

    }/*THEN*/

  }/*FOR*/

  if (!complete) {

    Free1D(missing);
    FreeGDT(sub_complete);

  }/*THEN*/

  FreeGDT(sub);

}/*RRD_GPERM*/

SEXP roundrobin_test(SEXP x, SEXP z, SEXP fixed, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP complete, SEXP debug) {

double *pvalue = NULL, a = NUM(alpha);
bool debugging = isTRUE(debug);
test_e test_type = test_to_enum(CHAR(STRING_ELT(test, 0)));
SEXP xx, zz, cc, which_fixed, result;

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, TRUE));
  PROTECT(zz = c_dataframe_column(data, z, FALSE, TRUE));
  /* match fixed variables. */
  PROTECT(which_fixed = match(fixed, z, 0));
  /* allocate the return value. */
  PROTECT(result = allocVector(REALSXP, length(z) - length(fixed)));
  setAttrib(result, R_NamesSymbol, string_setdiff(z, fixed));
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', length(result) * sizeof(double));

  /* extract the missing values indicators. */
  PROTECT(cc = subset_by_name(complete, 2, x, z));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    ddata dtx = ddata_from_SEXP(xx, 0), dtz = ddata_from_SEXP(zz, 0);

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 0, zz);
    meta_init_flags(&(dtz.m), 0, R_NilValue, which_fixed);

    rrd_discrete(dtx, dtz, test_type, pvalue, a, debugging);

    FreeDDT(dtx);
    FreeDDT(dtz);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    gdata dt = gdata_from_SEXP(zz, 1);

    meta_copy_names(&(dt.m), 1, zz);
    meta_init_flags(&(dt.m), 1, R_NilValue, which_fixed);
    dt.col[0] = REAL(VECTOR_ELT(xx, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));

    if (all_equal(cc, TRUESEXP)) {

      gdata_cache_means(&dt, 0);
      rrd_gaustests_complete(dt, test_type, pvalue, a, debugging);

    }/*THEN*/
    else {

      rrd_gaustests_with_missing(dt, test_type, pvalue, a, debugging);

    }/*ELSE*/

    FreeGDT(dt);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    cgdata dtx = cgdata_from_SEXP(xx, 0, 0), dtz = cgdata_from_SEXP(zz, 1, 1);

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 2, zz);
    meta_init_flags(&(dtz.m), 2, R_NilValue, which_fixed);

    if (all_equal(cc, TRUESEXP)) {

      rrd_micg_complete(dtx, dtz, test_type, pvalue, a, debugging);

    }/*THEN*/
    else {

      rrd_micg_with_missing(dtx, dtz, test_type, pvalue, a, debugging);

    }/*ELSE*/

    FreeCGDT(dtx);
    FreeCGDT(dtz);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    ddata dtx = ddata_from_SEXP(xx, 0), dtz = ddata_from_SEXP(zz, 0);

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 0, zz);
    meta_init_flags(&(dtz.m), 0, R_NilValue, which_fixed);

    rrd_dperm(dtx, dtz, test_type, pvalue, a, INT(B), IS_SMC(test_type) ? a : 1,
      debugging);

    FreeDDT(dtx);
    FreeDDT(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    gdata dt = gdata_from_SEXP(zz, 1);

    meta_copy_names(&(dt.m), 1, zz);
    meta_init_flags(&(dt.m), 1, R_NilValue, which_fixed);
    dt.col[0] = REAL(VECTOR_ELT(xx, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));

    rrd_gperm(dt, test_type, pvalue, a, INT(B), IS_SMC(test_type) ? a : 1,
      all_equal(cc, TRUESEXP), debugging);

    FreeGDT(dt);

  }/*THEN*/

  /* increment the test counter. */
  test_counter += length(zz) - length(fixed);

  UNPROTECT(5);

  return result;

}/*ROUNDROBIN_TEST*/
