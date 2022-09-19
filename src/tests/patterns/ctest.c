#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/sets.h"
#include "../../minimal/data.frame.h"
#include "../tests.h"
#include "../../core/moments.h"
#include "../../core/covariance.matrix.h"
#include "../../core/correlation.h"
#include "../../include/globals.h"
#include "../../core/data.table.h"
#include "../../math/linear.algebra.h"

/* parametric tests for discrete variables. */
double ct_discrete(ddata dtx, ddata dty, ddata dtz, double *pvalue, double *df,
    test_e test) {

int i = 0, llx = 0, lly = dty.nlvl[0], llz = 0;
int *xptr = NULL, *yptr = dty.col[0], *zptr = NULL;
double statistic = 0;

  zptr = Calloc1D(dtz.m.nobs, sizeof(int));
  c_fast_config(dtz.col, dtz.m.nobs, dtz.m.ncols, dtz.nlvl, zptr, &llz, 1);

  for (i = 0; i < dtx.m.ncols; i++) {

    xptr = dtx.col[i];
    llx = dtx.nlvl[i];

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, dtz.m.nobs,
                    df, test, (test == MI) || (test == MI_ADF));
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, dtz.m.nobs, df, TRUE);
      pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, dtz.m.nobs);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  Free1D(zptr);

  return statistic;

}/*CT_DISCRETE*/

/* parametric tests for Gaussian variables. */
double ct_gaustests_complete(gdata dtx, gdata dt, double *pvalue, double *df,
    test_e test) {

int i = 0, ntests = dtx.m.ncols;
double transform = 0, statistic = 0, lambda = 0;
covariance cov = { 0 }, basecov = { 0 };

  /* compute the degrees of freedom for correlation and mutual information. */
  *df = gaussian_cdf(test, dt.m.nobs, dt.m.ncols - 2);

  if (*df < 1) {

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
  c_covmat(dt.col, dt.mean, dt.m.nobs, dt.m.ncols, cov, 1);
  if (ntests > 1)
    copy_covariance(&cov, &basecov);

  /* compute the partial correlation and the test statistic. */
  for (i = 0; i < ntests; i++) {

    /* extract and plug-in the i-th variable. */
    dt.col[0] = dtx.col[i];
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
      statistic = cor_zf_trans(statistic, *df);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  FreeCOV(basecov);
  FreeCOV(cov);

  return statistic;

}/*CT_GAUSTESTS_COMPLETE*/

double ct_gaustests_with_missing(gdata dtx, gdata dt, double *pvalue,
    double *df, test_e test) {

int i = 0, ncomplete = 0;
double transform = 0, *mean = NULL, statistic = 0, lambda = 0;
bool *missing_yz = NULL, *missing_all = NULL;
covariance cov = { 0 };

  /* allocate the mean vector. */
  mean = Calloc1D(dt.m.ncols, sizeof(double));
  /* allocate the covariance matrix. */
  cov = new_covariance(dt.m.ncols, TRUE);

  /* allocate and initialize the missing values indicators. */
  missing_yz = Calloc1D(dt.m.nobs, sizeof(bool));
  missing_all = Calloc1D(dt.m.nobs, sizeof(bool));

  gdata_incomplete_cases(&dt, missing_yz, 1);

  /* compute the partial correlation and the test statistic. */
  for (i = 0; i < dtx.m.ncols; i++) {

    /* extract and plug in the i-th variable. */
    dt.col[0] = dtx.col[i];
    /* compute the covariance matrix. */
    c_covmat_with_missing(dt.col, dt.m.nobs, dt.m.ncols, missing_yz,
      missing_all, mean, cov.mat, &ncomplete);

    /* compute the degrees of freedom for correlation and mutual information. */
    *df = gaussian_cdf(test, ncomplete, dt.m.ncols - 2);

    if ((ncomplete == 0) || (*df < 1)) {

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
      statistic = cor_zf_trans(statistic, *df);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

  }/*FOR*/

  Free1D(mean);
  Free1D(missing_yz);
  Free1D(missing_all);
  FreeCOV(cov);

  return statistic;

}/*CT_GAUSTESTS_WITH_MISSING*/

/* conditional linear Gaussian mutual information test. */
double c_micg_chisq(cgdata dtx, cgdata dty, cgdata dtz, int *zptr, int llz,
    double *df, bool copy) {

double statistic = 0;

  if (dty.m.flag[0].discrete && dtx.m.flag[0].discrete) {

    if (dtz.ngcols - 1 > 0) {

      /* need to reverse conditioning to actually compute the test. */
      statistic = 2 * dtz.m.nobs * dtz.m.nobs *
                    c_cmicg_unroll(dtx.dcol[0], dtx.nlvl[0], dty.dcol[0],
                      dty.nlvl[0], zptr, llz, dtz.gcol + 1, dtz.ngcols - 1, df,
                      dtz.m.nobs);

    }/*THEN*/
    else {

      /* the test reverts back to a discrete mutual information test. */
      statistic = c_cchisqtest(dtx.dcol[0], dtx.nlvl[0],
                    dty.dcol[0], dty.nlvl[0], zptr, llz, dtz.m.nobs, df, MI, TRUE);

    }/*ELSE*/

  }/*THEN*/
  else if (dty.m.flag[0].gaussian && dtx.m.flag[0].gaussian) {

    if (copy)
      memcpy(dtz.gcol[0], dtx.gcol[0], dtz.m.nobs * sizeof(double));
    else
     dtz.gcol[0] = dtx.gcol[0];

    statistic = 2 * dtz.m.nobs * c_cmicg(dty.gcol[0], dtz.gcol, dtz.ngcols, NULL,
                      0, zptr, llz, dtz.nlvl, dtz.m.nobs, df);

  }/*THEN*/
  else if (dty.m.flag[0].discrete && dtx.m.flag[0].gaussian) {

    if (copy)
      memcpy(dtz.dcol[0], dty.dcol[0], dtz.m.nobs * sizeof(int));
    else
      dtz.dcol[0] = dty.dcol[0];

    dtz.nlvl[0] = dty.nlvl[0];
    statistic = 2 * dtz.m.nobs * c_cmicg(dtx.gcol[0], dtz.gcol + 1,
                      dtz.ngcols - 1, dtz.dcol, dtz.ndcols, zptr, llz, dtz.nlvl,
                      dtz.m.nobs, df);

  }/*THEN*/
  else if (dty.m.flag[0].gaussian && dtx.m.flag[0].discrete) {

    if (copy)
      memcpy(dtz.dcol[0], dtx.dcol[0], dtz.m.nobs * sizeof(int));
    else
      dtz.dcol[0] = dtx.dcol[0];

    dtz.nlvl[0] = dtx.nlvl[0];
    statistic = 2 * dtz.m.nobs * c_cmicg(dty.gcol[0], dtz.gcol + 1, dtz.ngcols - 1,
                      dtz.dcol, dtz.ndcols, zptr, llz, dtz.nlvl,
                      dtz.m.nobs, df);

  }/*ELSE*/

  return statistic;

}/*C_MICG_CHISQ*/

double ct_micg_complete(cgdata dtx, cgdata dty, cgdata dtz, double *pvalue,
    double *df) {

int i = 0, *zptr = 0, llz = 0;
double statistic = 0;
cgdata dtx_cur = { 0 };

  /* allocate the data table that will be filled with each variable in turn. */
  dtx_cur = empty_cgdata(dtx.m.nobs, 1, 1);

  /* allocate vector for the configurations of the discrete parents. */
  if (dtz.ndcols - 1 > 0) {

    zptr = Calloc1D(dtz.m.nobs, sizeof(int));
    c_fast_config(dtz.dcol + 1, dtz.m.nobs, dtz.ndcols - 1, dtz.nlvl + 1,
      zptr, &llz, 1);

  }/*THEN*/

  for (i = 0; i < dtx.m.ncols; i++) {

    /* extract the variable to test... */
    cgdata_subset_columns(&dtx, &dtx_cur, &i, 1);
    /* ... compute the test statistic... */
    statistic = c_micg_chisq(dtx_cur, dty, dtz, zptr, llz, df, FALSE);
    /* ... and save the p-value. */
    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);

  }/*FOR*/

  Free1D(zptr);
  FreeCGDT(dtx_cur);

  return statistic;

}/*CT_MICG_COMPLETE*/

double ct_micg_with_missing(cgdata dtx, cgdata dty, cgdata dt, double *pvalue,
    double *df) {

int i = 0, j = 0, *zptr = 0, llz = 0;
double statistic = 0;
bool *missing_x = NULL, *missing_yz = NULL, *missing_all = NULL;
cgdata dt_complete = { 0 }, dtx_cur = { 0 }, dtx_complete = { 0 };
cgdata dty_complete = { 0 };

  /* allocate and initialize a data table for the complete observations. */
  dt_complete = new_cgdata(dt.m.nobs, dt.ndcols, dt.ngcols);
  dtx_cur = empty_cgdata(dtx.m.nobs, 1, 1);
  dtx_complete = new_cgdata(dtx.m.nobs, dtx.ndcols, dtx.ngcols);
  dty_complete = new_cgdata(dty.m.nobs, dty.ndcols, dty.ngcols);

  /* compute the missing values indicators. */
  missing_x = Calloc1D(dt.m.nobs, sizeof(bool));
  missing_yz = Calloc1D(dt.m.nobs, sizeof(bool));
  missing_all = Calloc1D(dt.m.nobs, sizeof(bool));

  cgdata_incomplete_cases(&dty, missing_yz, 0, 0);
  cgdata_incomplete_cases(&dt, missing_yz, 1, 1);

  for (i = 0; i < dtx.m.ncols; i++) {

    /* extract the variable to test ... */
    cgdata_subset_columns(&dtx, &dtx_cur, &i, 1);
    /* ... find out which observations are missing ... */
    cgdata_incomplete_cases(&dtx_cur, missing_x, 0, 0);
    /* ... update the global missingness indicators ... */
    for (j = 0; j < dt.m.nobs; j++)
      missing_all[j] = missing_x[j] || missing_yz[j];
    /* ... and subset the complete data. */
    cgdata_subsample_by_logical(&dtx_cur, &dtx_complete, missing_all, 0, 0);
    cgdata_subsample_by_logical(&dt, &dt_complete, missing_all, 1, 1);
    cgdata_subsample_by_logical(&dty, &dty_complete, missing_all, 0, 0);

    /* assume independence and return if there are no complete observations. */
    if (dt_complete.m.nobs == 0) {

      statistic = 0;
      pvalue[i] = 1;
      continue;

    }/*THEN*/

    if (dt.ndcols - 1 > 0) {

     zptr = Calloc1D(dt.m.nobs, sizeof(int));
      c_fast_config(dt_complete.dcol + 1, dt_complete.m.nobs,
        dt_complete.ndcols - 1, dt_complete.nlvl + 1, zptr, &llz, 1);

    }/*THEN*/

    statistic = c_micg_chisq(dtx_complete, dty_complete, dt_complete,
                  zptr, llz, df, TRUE);

    pvalue[i] = pchisq(statistic, *df, FALSE, FALSE);
    Free1D(zptr);

  }/*FOR*/

  Free1D(missing_x);
  Free1D(missing_yz);
  Free1D(missing_all);
  FreeCGDT(dt_complete);
  FreeCGDT(dtx_cur);
  FreeCGDT(dtx_complete);
  FreeCGDT(dty_complete);

  return statistic;

}/*CT_MICG_WITH_MISSING*/

/* discrete permutation tests. */
double ct_dperm(ddata dtx, ddata dty, ddata dtz, double *pvalue, double *df,
    test_e type, int B, double a) {

int i = 0, *zptr = NULL, llz = 0;
double statistic = 0;

  zptr = Calloc1D(dtz.m.nobs, sizeof(int));
  c_fast_config(dtz.col, dtz.m.nobs, dtz.m.ncols, dtz.nlvl, zptr, &llz, 1);

  for (i = 0; i < dtx.m.ncols; i++) {

    c_cmcarlo(dtx.col[i], dtx.nlvl[i], dty.col[0], dty.nlvl[0], zptr, llz,
      dtx.m.nobs, B, &statistic, pvalue + i, a, type, df);

  }/*FOR*/

  Free1D(zptr);

  return statistic;

}/*CT_DPERM*/

/* continuous permutation tests. */
double ct_gperm(gdata dtx, gdata dt, double *pvalue, double *df, test_e type,
    int B, double a, bool complete) {

int i = 0, j = 0, k = 0, nc = 0;
double *yptr = dt.col[1], **complete_column = NULL, statistic = 0;
bool *missing_yz = NULL;
gdata dt_complete = { 0 };

  if (!complete) {

    /* allocate and initialize the missing values indicators. */
    missing_yz = Calloc1D(dt.m.nobs, sizeof(bool));
    gdata_incomplete_cases(&dt, missing_yz, 1);

    dt_complete = new_gdata(dt.m.nobs, dt.m.ncols);
    complete_column = dt_complete.col;

  }/*THEN*/
  else {

    complete_column = dt.col;

  }/*ELSE*/

  for (i = 0; i < dtx.m.ncols; i++) {

    /* swap the first column and restore the second, which is that undergoing
     * permutation (backward compatibility from set random seed). */
    dt.col[0] = dtx.col[i];
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

    c_gauss_cmcarlo(complete_column, dt.m.ncols, nc, 0, 1, B, &statistic,
      pvalue + i, a, type);

  }/*FOR*/

  if (!complete) {

    Free1D(missing_yz);
    FreeGDT(dt_complete);

  }/*THEN*/

  return statistic;

}/*CT_GPERM*/

