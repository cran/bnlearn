#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/data.table.h"
#include "include/matrix.h"

static void rrd_disc_message(meta m, SEXP x, int offset, const char *y,
    double p, double alpha) {

  Rprintf("    > node %s is %s %s given ", CHAR(STRING_ELT(x, 0)),
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
static void rrd_discrete(SEXP xx, SEXP zz, SEXP fixed, SEXP x, test_e test,
    double *pvalue, double alpha, int debuglevel) {

int i = 0, cur = 0, valid = 0, llx = NLEVELS(xx), lly = 0, llz = 0;
int *xptr = INTEGER(xx), *yptr = NULL, *zptr = NULL;
double statistic = 0, df = 0;
ddata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = ddata_from_SEXP(zz, 0);
  meta_copy_names(&(dt.m), 0, zz);
  meta_init_flags(&(dt.m), 0, R_NilValue, fixed);
  valid = dt.m.ncols;
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dt.m.nobs, sizeof(int));

  for (i = 0; i < dt.m.ncols; i++) {

    /* nothing to do if there is just one variable left. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* extract the variable to test. */
    yptr = dt.col[i];
    lly = dt.nlvl[i];
    /* drop the variable to test from the conditioning set. */
    dt.m.flag[i].drop = TRUE;
    ddata_drop_flagged(&dt, &sub);
    /* construct the configurations of the conditioning variables. */
    c_fast_config(sub.col, sub.m.nobs, sub.m.ncols, sub.nlvl, zptr, &llz, 1);

    if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

      /* mutual information and Pearson's X^2 asymptotic tests. */
      statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs,
                    &df, test, (test == MI) || (test == MI_ADF));
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == MI_SH) {

      /* shrinkage mutual information test. */
      statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, &df, TRUE);
      pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    }/*THEN*/
    else if (test == JT) {

      /* Jonckheere-Terpstra test. */
      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs);
      pvalue[cur] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    if (debuglevel > 0)
      rrd_disc_message(sub.m, x, 0, dt.m.names[i], pvalue[cur], alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dt.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(dt, FALSE);
  FreeDDT(sub, FALSE);

}/*RRD_DISCRETE*/

/* parametric tests for gaussian variables. */
static void rrd_gaustests(SEXP xx, SEXP zz, SEXP fixed, SEXP x, test_e test,
    double *pvalue, double alpha, int debuglevel) {

int i = 1, cur = 0, df = 0, valid = 0, t = 0, run_svd = TRUE;
double statistic = 0, lambda = 0;
gdata dt = { 0 }, sub = { 0 };
covariance cov = { 0 }, backup = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 1);
  meta_copy_names(&(dt.m), 1, zz);
  meta_init_flags(&(dt.m), 1, R_NilValue, fixed);
  valid = dt.m.ncols - 1;
  dt.col[0] = REAL(xx);
  dt.m.names[0] = CHAR(STRING_ELT(x, 0));
  gdata_cache_means(&dt, 0);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.m.names = Calloc1D(dt.m.ncols, sizeof(char *));
  sub.m.flag = Calloc1D(dt.m.ncols, sizeof(flags));
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
    if (test == COR)
      df = sub.m.nobs - sub.m.ncols;
    else if ((test == MI_G) || (test == MI_G_SH))
      df = 1;

    if (((test == COR) && (df < 1)) ||
        ((test == ZF) && (sub.m.nobs - sub.m.ncols < 2))) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      pvalue[cur++] = 1;

      if (debuglevel > 0)
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
      statistic = cor_zf_trans(statistic, (double)(sub.m.nobs - sub.m.ncols));
      pvalue[cur] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*THEN*/

    if (debuglevel > 0)
      rrd_gauss_message(sub, t, pvalue[cur], alpha);

    /* restore the covariance matrix before subsetting it. */
    copy_covariance(&backup, &cov);

    /* remove ifrom the covariance matrix and make sure to recompute SVD. */
    if (pvalue[cur++] > alpha) {

      valid--;
      dt.m.flag[i].drop = TRUE;
      covariance_drop_column(&cov, t);
      run_svd = TRUE;

    }/*THEN*/
    else {

      run_svd = FALSE;

    }/*ELSE*/

  }/*WHILE*/

  FreeGDT(dt, FALSE);
  FreeGDT(sub, FALSE);
  FreeCOV(backup);
  FreeCOV(cov);

}/*RRD_GAUSTESTS*/

/* conditional linear Gaussian test. */
static void rrd_micg(SEXP xx, SEXP zz, SEXP fixed, SEXP x, test_e test,
    double *pvalue, double alpha, int debuglevel) {

int i = 0;
int cur = 0, valid = 0, xtype = TYPEOF(xx), llx = 0, lly = 0, llz = 0;
int *configurations = NULL, *zptr = NULL;
double statistic = 0, df = 0;
void *xptr = NULL, *yptr = NULL;
cgdata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = cgdata_from_SEXP(zz, 1, 1);
  meta_copy_names(&(dt.m), 2, zz);
  meta_init_flags(&(dt.m), 2, R_NilValue, fixed);
  valid = dt.m.ncols - 2;
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_cgdata(dt.m.nobs, dt.ndcols, dt.ngcols);

  /* dereference xx, whatever the type. */
  if (xtype == INTSXP) {

    xptr = INTEGER(xx);
    llx = NLEVELS(xx);

  }/*THEN*/
  else {

    xptr = REAL(xx);

  }/*ELSE*/

  configurations = Calloc1D(dt.m.nobs, sizeof(int));

  for (i = 2; i < dt.m.ncols; i++) {

    /* nothing to do if there are no more valid nodes.*/
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* drop the variable to test from the conditioning set. */
    dt.m.flag[i].drop = TRUE;
    cgdata_drop_flagged(&dt, &sub);

    /* extract the variable to test. */
    if (dt.m.flag[i].discrete) {

      yptr = dt.dcol[dt.map[i]];
      lly = dt.nlvl[dt.map[i]];

    }/*THEN*/
    else {

      yptr = dt.gcol[dt.map[i]];

    }/*ELSE*/

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

    if ((xtype == INTSXP) && dt.m.flag[i].discrete) {

      /* check whether the conditioning set is valid. */
      if (sub.ngcols > 1) {

        /* need to reverse conditioning to actually compute the test. */
        statistic = 2 * sub.m.nobs * sub.m.nobs *
                      c_cmicg_unroll(xptr, llx, yptr, lly, zptr, llz,
                                 sub.gcol + 1, sub.ngcols - 1, &df, sub.m.nobs);

      }/*THEN*/
      else {

        /* if both nodes are discrete, the test reverts back to a discrete
         * mutual information test. */
        statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs,
                      &df, MI, TRUE);

      }/*ELSE*/

    }/*THEN*/
    else if ((xtype == REALSXP) && dt.m.flag[i].gaussian) {

      sub.gcol[0] = xptr;
      statistic = 2 * sub.m.nobs * c_cmicg(yptr, sub.gcol, sub.ngcols, NULL, 0,
                                    zptr, llz, 0, sub.m.nobs, &df);

    }/*THEN*/
    else if ((xtype == REALSXP) && dt.m.flag[i].discrete) {

      sub.dcol[0] = yptr;
      sub.nlvl[0] = lly;
      statistic = 2 * sub.m.nobs *
                    c_cmicg(xptr, sub.gcol + 1, sub.ngcols - 1, sub.dcol,
                             sub.ndcols, zptr, llz, sub.nlvl, sub.m.nobs, &df);

    }/*THEN*/
    else if ((xtype == INTSXP) && dt.m.flag[i].gaussian) {

      sub.dcol[0] = xptr;
      sub.nlvl[0] = llx;
      statistic = 2 * sub.m.nobs *
                    c_cmicg(yptr, sub.gcol + 1, sub.ngcols - 1, sub.dcol,
                             sub.ndcols, zptr, llz, sub.nlvl, sub.m.nobs, &df);

    }/*THEN*/

    pvalue[cur] = pchisq(statistic, df, FALSE, FALSE);

    if (debuglevel > 0)
      rrd_disc_message(sub.m, x, 2, dt.m.names[i], pvalue[cur], alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dt.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(configurations);
  FreeCGDT(dt, FALSE);
  FreeCGDT(sub, FALSE);

}/*RRD_MICG*/

/* discrete permutation tests. */
static void rrd_dperm(SEXP xx, SEXP zz, SEXP fixed, SEXP x, test_e test,
    double *pvalue, double alpha, int nperms, double threshold, int debuglevel) {

int i = 0, cur = 0, valid = 0, llx = NLEVELS(xx), lly = 0, llz = 0;
int *xptr = INTEGER(xx), *yptr = NULL, *zptr = NULL;
double statistic = 0, df = 0;
ddata dt = {0}, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = ddata_from_SEXP(zz, 0);
  meta_copy_names(&(dt.m), 0, zz);
  meta_init_flags(&(dt.m), 0, R_NilValue, fixed);
  valid = dt.m.ncols;
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dt.m.nobs, sizeof(int));

  for (i = 0; i < dt.m.ncols; i++) {

    /* nothing to do if there is just one variable. */
    if (valid < 2)
      break;
    /* do not test fixed variables. */
    if (dt.m.flag[i].fixed)
      continue;

    /* extract the variable to test. */
    yptr = dt.col[i];
    lly = dt.nlvl[i];
    /* drop the variable to test from the conditioning set. */
    dt.m.flag[i].drop = TRUE;
    ddata_drop_flagged(&dt, &sub);
    /* construct the configurations of the conditioning variables. */
    c_fast_config(sub.col, sub.m.nobs, sub.m.ncols, sub.nlvl, zptr, &llz, 1);

    c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, nperms, &statistic,
      pvalue + cur, threshold, test, &df);

    if (debuglevel > 0)
      rrd_disc_message(sub.m, x, 0, dt.m.names[i], pvalue[cur], alpha);

    /* do not drop the variable if the tests is significant. */
    if (pvalue[cur++] > alpha)
      valid--;
    else
      dt.m.flag[i].drop = FALSE;

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(dt, FALSE);
  FreeDDT(sub, FALSE);

}/*RRD_DPERM*/

/* continuous permutation tests. */
static void rrd_gperm(SEXP xx, SEXP zz, SEXP fixed, SEXP x, test_e test,
    double *pvalue, double alpha, int nperms, double threshold, int debuglevel) {

int i = 0, cur = 0, valid = 0, t = 0;
double statistic = 0;
gdata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 1);
  meta_copy_names(&(dt.m), 1, zz);
  meta_init_flags(&(dt.m), 1, R_NilValue, fixed);
  valid = dt.m.ncols - 1;
  dt.col[0] = REAL(xx);
  dt.m.names[0] = CHAR(STRING_ELT(x, 0));

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.m.names = Calloc1D(dt.m.ncols, sizeof(char *));
  sub.m.flag = Calloc1D(dt.m.ncols, sizeof(flags));

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

    c_gauss_cmcarlo(sub.col, sub.m.ncols, sub.m.nobs, 0, t, nperms,
      &statistic, pvalue + cur, threshold, test);

    if (debuglevel > 0)
      rrd_gauss_message(sub, t, pvalue[cur], alpha);

    /* remove from the set of valid candidates for the conditioning set. */
    if (pvalue[cur++] > alpha) {

      valid--;
      dt.m.flag[i].drop = TRUE;

    }/*THEN*/

  }/*FOR*/

  FreeGDT(dt, FALSE);
  FreeGDT(sub, FALSE);

}/*RRD_GPERM*/

SEXP roundrobin_test(SEXP x, SEXP z, SEXP fixed, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP complete, SEXP debug) {

int debuglevel = isTRUE(debug);
double *pvalue = NULL, a = NUM(alpha);
test_e test_type = test_label(CHAR(STRING_ELT(test, 0)));
SEXP xx, zz, which_fixed, result;

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, z, FALSE, TRUE));
  /* match fixed variables. */
  PROTECT(which_fixed = match(fixed, z, 0));
  /* allocate the return value. */
  PROTECT(result = allocVector(REALSXP, length(z) - length(fixed)));
  setAttrib(result, R_NamesSymbol, string_setdiff(z, fixed));
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', length(result) * sizeof(double));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    rrd_discrete(xx, zz, which_fixed, x, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    rrd_gaustests(xx, zz, which_fixed, x, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    rrd_micg(xx, zz, which_fixed, x, test_type, pvalue, a, debuglevel);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    rrd_dperm(xx, zz, which_fixed, x, test_type, pvalue, a, INT(B),
      IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    rrd_gperm(xx, zz, which_fixed, x, test_type, pvalue, a, INT(B),
      IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/

  /* increment the test counter. */
  test_counter += length(zz);

  UNPROTECT(4);

  return result;

}/*ROUNDROBIN_TEST*/
