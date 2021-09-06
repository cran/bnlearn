#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/data.table.h"

static void update_pvalue_range(double pvalue, double *min, double *max) {

  *min = pvalue < *min ? pvalue : *min;
  *max = pvalue > *max ? pvalue : *max;

}/*UPDATE_PVALUES*/

static SEXP ast_prepare_retval(double pvalue, double min_pvalue,
  double max_pvalue, double alpha, const char **nodes, int nnodes) {

int i = 0;
SEXP retval, dsep_set;

  PROTECT(retval = mkRealVec(3, pvalue, min_pvalue, max_pvalue));
  setAttrib(retval, R_NamesSymbol,
    mkStringVec(3, "p.value", "min.p.value", "max.p.value"));

  if (pvalue > alpha) {

    PROTECT(dsep_set = allocVector(STRSXP, nnodes));
    for (i = 0; i < nnodes; i++)
      SET_STRING_ELT(dsep_set, i, mkChar(nodes[i]));
    setAttrib(retval, BN_DsepsetSymbol, dsep_set);
    UNPROTECT(1);

  }/*THEN*/
  else {

    setAttrib(retval, BN_DsepsetSymbol, R_NilValue);

  }/*ELSE*/

  UNPROTECT(1);

  return retval;

}/*AST_PREPARE_RETVAL*/

/* parametric tests for discrete variables. */
static SEXP ast_discrete(ddata dtx, ddata dty, ddata dtz, int nf, int minsize,
    int maxsize, test_e test, double a, bool debugging) {

int *xptr = dtx.col[0], *yptr = dty.col[0], *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = dtx.nlvl[0], lly = dty.nlvl[0], llz = 0;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;
SEXP retval;
ddata sub = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dtz.m.nobs, dtz.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dtz.m.nobs, sizeof(int));

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset. */
    subset = Calloc1D(cursize + nf, sizeof(int));
    /* initialize the first subset (assumes fixed variables are the first nf). */
    first_subset(subset + nf, cursize, nf);
    for (i = 0; i < nf; i++)
      subset[i] = i;

    /* iterate over subsets. */
    do {

      /* prepare the current subset. */
      ddata_subset_columns(&dtz, &sub, subset, cursize + nf);
      /* construct the parents' configurations. */
      c_fast_config(sub.col, sub.m.nobs, cursize + nf, sub.nlvl, zptr, &llz, 1);

      if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

        /* mutual information and Pearson's X^2 asymptotic tests. */
        statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, &df,
                      test, (test == MI) || (test == MI_ADF));
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == MI_SH) {

        /* shrinkage mutual information test. */
        statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, &df, TRUE);
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == JT) {

        /* Jonckheere-Terpstra test. */
        statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs);
        pvalue = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ", dtx.m.names[0],
          (pvalue > a) ? "independent from" : "dependent on", dty.m.names[0]);
        for (i = 0; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names, sub.m.ncols));

        Free1D(subset);
        Free1D(zptr);
        FreeDDT(sub);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf, cursize, dtz.m.ncols - nf, nf));

    Free1D(subset);

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_DISCRETE*/

/* parametric tests for Gaussian variables (for complete data). */
static SEXP ast_gaustests_complete(gdata dt, int nf, int minsize, int maxsize,
    double a, bool debugging, test_e test) {

int i = 0, cursize = 0, *subset = NULL, df = 0;
double statistic = 0, lambda = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
SEXP retval;
gdata sub = { 0 };
covariance cov = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* compute the degrees of freedom for correlation and mutual information. */
    df = gaussian_cdf(test, dt.m.nobs, cursize + nf);

    if (df < 1) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      pvalue = min_pvalue = max_pvalue = 1;

      if (debugging) {

        Rprintf("    > node %s is independent from %s given any conditioning set of size %d ",
          dt.m.names[0], dt.m.names[1], cursize + nf);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      /* the conditioning set is the first possible conditioning set, which
       * comprises the first columns of the data table. */
      PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                         a, dt.m.names + 2, cursize + nf));

      FreeGDT(sub);

      UNPROTECT(1);
      return retval;

    }/*THEN*/

    /* allocate and initialize the subset indexes array. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* allocate the covariance matrix and the U, D, V matrix. */
    cov = new_covariance(cursize + nf + 2, TRUE);
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    do {

      /* prepare the current subset. */
      gdata_subset_columns(&dt, &sub, subset, cursize + nf + 2);
      /* compute the covariance matrix. */
      c_covmat(sub.col, sub.mean, sub.m.nobs, sub.m.ncols, cov, 0);

      if (test == COR) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = cor_t_trans(statistic, (double)df);
        pvalue = 2 * pt(fabs(statistic), df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == MI_G) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == MI_G_SH) {

        lambda = covmat_lambda(sub.col, sub.mean, cov, sub.m.nobs, NULL,
                   sub.m.nobs);
        covmat_shrink(cov, lambda);
        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == ZF) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = cor_zf_trans(statistic, df);
        pvalue = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ",
          sub.m.names[0], (pvalue > a) ? "independent from" : "dependent on",
          sub.m.names[1]);
        for (i = 2; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names + 2, sub.m.ncols - 2));

        Free1D(subset);
        FreeCOV(cov);
        FreeGDT(sub);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    FreeCOV(cov);
    Free1D(subset);

  }/*FOR*/

  FreeGDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_GAUSTESTS_COMPLETE*/

/* parametric tests for Gaussian variables (for incomplete data). */
static SEXP ast_gaustests_with_missing(gdata dt, int nf, int minsize,
    int maxsize, double a, bool debugging, test_e test) {

int i = 0, cursize = 0, *subset = NULL, df = 0, ncomplete = 0;
double statistic = 0, lambda = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
double *mean = NULL;
bool *missing_xy = NULL, *missing_all = NULL;
SEXP retval;
gdata sub = { 0 };
covariance cov = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));
  /* allocate the missingness indicators. */
  missing_xy = Calloc1D(dt.m.nobs, sizeof(bool));
  gdata_incomplete_cases_range(&dt, missing_xy, 0, 1);

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate a vector to store column means (from complete observations). */
    mean = Calloc1D(cursize + nf + 2, sizeof(double));
    /* allocate missingness indicators. */
    missing_all = Calloc1D(sub.m.nobs, sizeof(bool));
    /* allocate and initialize the subset indexes array. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* allocate the covariance matrix and the U, D, V matrix. */
    cov = new_covariance(cursize + nf + 2, TRUE);
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    do {

      /* prepare the current subset. */
      gdata_subset_columns(&dt, &sub, subset, cursize + nf + 2);
      /* compute the covariance matrix. */
      c_covmat_with_missing(sub.col, sub.m.nobs, sub.m.ncols, missing_xy,
        missing_all, mean, cov.mat, &ncomplete);

      /* compute the degrees of freedom for correlation and mutual information. */
      df = gaussian_cdf(test, ncomplete, cursize + nf);

      if ((ncomplete == 0) || (df < 1)) {

        /* if there are not enough degrees of freedom, return independence. */
        warning("trying to do a conditional independence test with zero degrees of freedom.");

        pvalue = min_pvalue = max_pvalue = 1;

        if (debugging) {

          Rprintf("    > node %s is independent from %s given any conditioning set of size %d ",
            dt.m.names[0], dt.m.names[1], cursize + nf);
          Rprintf("(p-value: %g).\n", pvalue);

        }/*THEN*/

        /* the conditioning set is the first possible conditioning set, which
         * comprises the first columns of the data table. */
        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, dt.m.names + 2, cursize + nf));

        FreeGDT(sub);
        Free1D(missing_xy);
        Free1D(missing_all);
        Free1D(subset);
        Free1D(mean);
        FreeCOV(cov);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

      if (test == COR) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = cor_t_trans(statistic, (double)df);
        pvalue = 2 * pt(fabs(statistic), df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == MI_G) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * ncomplete * cor_mi_trans(statistic);
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == MI_G_SH) {

        lambda = covmat_lambda(sub.col, mean, cov, sub.m.nobs, missing_all,
                   ncomplete);
        covmat_shrink(cov, lambda);
        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * ncomplete * cor_mi_trans(statistic);
        pvalue = pchisq(statistic, df, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/
      else if (test == ZF) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = cor_zf_trans(statistic, df);
        pvalue = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);
        update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      }/*THEN*/

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ",
          sub.m.names[0], (pvalue > a) ? "independent from" : "dependent on",
          sub.m.names[1]);
        for (i = 2; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (FALSE &&  pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names + 2, sub.m.ncols - 2));

        Free1D(subset);
        FreeCOV(cov);
        FreeGDT(sub);
        Free1D(mean);
        Free1D(missing_xy);
        Free1D(missing_all);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    FreeCOV(cov);
    Free1D(subset);
    Free1D(mean);
    Free1D(missing_all);

  }/*FOR*/

  Free1D(missing_xy);
  FreeGDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_GAUSTESTS_WITH_MISSING*/

/* conditional linear Gaussian test (for complete data). */
static SEXP ast_micg_complete(cgdata dtx, cgdata dty, cgdata dtz, int nf,
    int minsize, int maxsize, double a, bool debugging) {

int i = 0, *subset = NULL, cursize = 0;
int *zptr = NULL, llx = 0, lly = 0, llz = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
void *xptr = 0, *yptr = 0;
SEXP retval;
cgdata sub = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);

  /* if both variables are continuous and all conditioning variables are
   * continuous, the test reverts back to a Gaussian mutual information test. */
  if (dtx.m.flag[0].discrete && dty.m.flag[0].discrete) {

    xptr = dtx.dcol[0];
    llx = dtx.nlvl[0];
    yptr = dty.dcol[0];
    lly = dty.nlvl[0];

  }/*THEN*/
  else if (dtx.m.flag[0].gaussian && dty.m.flag[0].gaussian) {

    xptr = dtx.gcol[0];
    yptr = dty.gcol[0];

  }/*THEN*/
  else if (dtx.m.flag[0].gaussian && dty.m.flag[0].discrete) {

    xptr = dtx.gcol[0];
    yptr = dty.dcol[0];
    lly = dty.nlvl[0];

  }/*THEN*/
  else if (dtx.m.flag[0].discrete && dty.m.flag[0].gaussian) {

    yptr = dtx.dcol[0];
    lly = dtx.nlvl[0];
    xptr = dty.gcol[0];

  }/*THEN*/

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    do {

      /* prepare the current subset. */
      cgdata_subset_columns(&dtz, &sub, subset, cursize + nf + 2);

      /* if there are discrete conditioning variables, compute their
       * configurations. */
      if (sub.ndcols - 1 > 0) {

        zptr = Calloc1D(sub.m.nobs, sizeof(int));
        c_fast_config(sub.dcol + 1, sub.m.nobs, sub.ndcols - 1, sub.nlvl + 1,
          zptr, &llz, 1);

      }/*THEN*/
      else {

        zptr = NULL;
        llz = 0;

      }/*ELSE*/

      if (dtx.m.flag[0].discrete && dty.m.flag[0].discrete) {

        /* check whether the conditioning set is valid. */
        if (sub.ngcols - 1 > 0) {

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
      else if (dtx.m.flag[0].gaussian && dty.m.flag[0].gaussian) {

        sub.gcol[0] = xptr;
        statistic = 2 * sub.m.nobs *
                      c_cmicg(yptr, sub.gcol, sub.ngcols, NULL, 0, zptr, llz,
                        sub.nlvl, sub.m.nobs, &df);

      }/*THEN*/
      else { /* one variable is discrete, the other is continuous. */

        sub.dcol[0] = yptr;
        sub.nlvl[0] = lly;
        statistic = 2 * sub.m.nobs *
          c_cmicg(xptr, sub.gcol + 1, sub.ngcols - 1, sub.dcol, sub.ndcols,
            zptr, llz, sub.nlvl, sub.m.nobs, &df);

      }/*ELSE*/

      Free1D(zptr);

      pvalue = pchisq(statistic, df, FALSE, FALSE);
      update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ", dtx.m.names[0],
          (pvalue > a) ? "independent from" : "dependent on", dty.m.names[0]);
        for (i = 2; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names + 2, sub.m.ncols - 2));

        Free1D(subset);
        FreeCGDT(sub);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dtz.m.ncols - nf - 2, nf + 2));

    Free1D(subset);

  }/*FOR*/

  FreeCGDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_MICG_COMPLETE*/

/* conditional linear Gaussian test (for incomplete data). */
static SEXP ast_micg_with_missing(cgdata dtx, cgdata dty, cgdata dtz, int nf,
    int minsize, int maxsize, double a, bool debugging) {

int i = 0, *subset = NULL, cursize = 0;
int *zptr = NULL, llz = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
bool *missing_xy = NULL, *missing_all = NULL;
SEXP retval;
cgdata sub = { 0 }, sub_complete = { 0 };
cgdata dtx_complete = { 0 }, dty_complete = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);
  sub_complete = new_cgdata(dtz.m.nobs, dtz.ndcols, dtz.ngcols);
  sub_complete.m.names = Calloc1D(sub_complete.m.ncols, sizeof(char *));
  dtx_complete = new_cgdata(dtx.m.nobs, dtx.ndcols, dtx.ngcols);
  dty_complete = new_cgdata(dty.m.nobs, dty.ndcols, dty.ngcols);

  /* allocate missingness indicators. */
  missing_xy = Calloc1D(dtz.m.nobs, sizeof(bool));
  missing_all = Calloc1D(dtz.m.nobs, sizeof(bool));

  cgdata_incomplete_cases(&dtx, missing_xy, 0, 0);
  cgdata_incomplete_cases(&dty, missing_xy, 0, 0);

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    do {

      /* extract the variables in the conditioning set ... */
      cgdata_subset_columns(&dtz, &sub, subset, cursize + nf + 2);
      /* ... find out which observations are missing ... */
      memcpy(missing_all, missing_xy, dtz.m.nobs * sizeof(bool));
      cgdata_incomplete_cases(&sub, missing_all, 1, 1);
      /*  ... and subset the complete data. */
      cgdata_subsample_by_logical(&sub, &sub_complete, missing_all, 1, 1);
      cgdata_subsample_by_logical(&dtx, &dtx_complete, missing_all, 0, 0);
      cgdata_subsample_by_logical(&dty, &dty_complete, missing_all, 0, 0);

      /* assume independence and return if there are no complete observations. */
      if (sub_complete.m.nobs <= 1) {

        pvalue = min_pvalue = max_pvalue = 1;
        goto exit;

      }/*THEN*/

      /* if there are discrete conditioning variables, compute their
       * configurations. */
      if (sub_complete.ndcols - 1 > 0) {

        zptr = Calloc1D(sub_complete.m.nobs, sizeof(int));
        c_fast_config(sub_complete.dcol + 1, sub_complete.m.nobs,
          sub_complete.ndcols - 1, sub_complete.nlvl + 1, zptr, &llz, 1);

      }/*THEN*/
      else {

        zptr = NULL;
        llz = 0;

      }/*ELSE*/

      if (dtx.m.flag[0].discrete && dty.m.flag[0].discrete) {

        /* check whether the conditioning set is valid. */
        if (sub_complete.ngcols - 1 > 0) {

          /* need to reverse conditioning to actually compute the test. */
          statistic = 2 * sub_complete.m.nobs * sub_complete.m.nobs *
                        c_cmicg_unroll(dtx_complete.dcol[0],
                          dtx_complete.nlvl[0], dty_complete.dcol[0],
                          dty_complete.nlvl[0], zptr, llz, sub_complete.gcol + 1,
                          sub_complete.ngcols - 1, &df, sub_complete.m.nobs);

        }/*THEN*/
        else {

          /* if both nodes are discrete, the test reverts back to a discrete
           * mutual information test. */
          statistic = c_cchisqtest(dtx_complete.dcol[0], dtx_complete.nlvl[0],
                        dty_complete.dcol[0], dty_complete.nlvl[0], zptr, llz,
                        sub_complete.m.nobs, &df, MI, TRUE);

        }/*ELSE*/

      }/*THEN*/
      else if (dtx.m.flag[0].gaussian && dty.m.flag[0].gaussian) {

        memcpy(sub_complete.gcol[0], dtx_complete.gcol[0],
          sub_complete.m.nobs * sizeof(double));
        statistic = 2 * sub_complete.m.nobs *
                      c_cmicg(dty_complete.gcol[0], sub_complete.gcol,
                        sub_complete.ngcols, NULL, 0, zptr, llz,
                        sub_complete.nlvl, sub_complete.m.nobs, &df);

      }/*THEN*/
      else if (dtx.m.flag[0].gaussian && dty.m.flag[0].discrete) {

        memcpy(sub_complete.dcol[0], dty_complete.dcol[0],
          sub_complete.m.nobs * sizeof(int));
        sub_complete.nlvl[0] = dty_complete.nlvl[0];
        statistic = 2 * sub_complete.m.nobs *
          c_cmicg(dtx_complete.gcol[0], sub_complete.gcol + 1,
            sub_complete.ngcols - 1, sub_complete.dcol, sub_complete.ndcols,
            zptr, llz, sub_complete.nlvl, sub_complete.m.nobs, &df);

      }/*THEN*/
      else if (dtx.m.flag[0].discrete && dty.m.flag[0].gaussian) {

        memcpy(sub_complete.dcol[0], dtx_complete.dcol[0],
          sub_complete.m.nobs * sizeof(int));
        sub_complete.nlvl[0] = dtx_complete.nlvl[0];
        statistic = 2 * sub_complete.m.nobs *
          c_cmicg(dty_complete.gcol[0], sub_complete.gcol + 1,
            sub_complete.ngcols - 1, sub_complete.dcol, sub_complete.ndcols,
            zptr, llz, sub_complete.nlvl, sub_complete.m.nobs, &df);

      }/*ELSE*/

      pvalue = pchisq(statistic, df, FALSE, FALSE);
      update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      /* increment the test counter. */
      test_counter++;

exit:

      Free1D(zptr);

      if (debugging) {

        Rprintf("    > node %s is %s %s given ", dtx.m.names[0],
          (pvalue > a) ? "independent from" : "dependent on", dty.m.names[0]);
        for (i = 2; i < sub_complete.m.ncols; i++)
          Rprintf("%s ", sub_complete.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub_complete.m.names + 2, sub_complete.m.ncols - 2));

        Free1D(subset);
        FreeCGDT(sub);
        FreeCGDT(sub_complete);
        FreeCGDT(dtx_complete);
        FreeCGDT(dty_complete);
        Free1D(missing_xy);
        Free1D(missing_all);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dtz.m.ncols - nf - 2, nf + 2));

    Free1D(subset);

  }/*FOR*/

  FreeCGDT(sub);
  FreeCGDT(sub_complete);
  FreeCGDT(dtx_complete);
  FreeCGDT(dty_complete);
  Free1D(missing_xy);
  Free1D(missing_all);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_MICG_WITH_MISSING*/

/* discrete permutation tests. */
static SEXP ast_dperm(ddata dtx, ddata dty, ddata dtz, int nf, int minsize,
    int maxsize, double a, test_e type, int nperms, double threshold,
    bool debugging) {

int *xptr = dtx.col[0], *yptr = dty.col[0], *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = dtx.nlvl[0], lly = dty.nlvl[0], llz = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
SEXP retval;
ddata sub = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dtz.m.nobs, dtz.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dtz.m.nobs, sizeof(int));

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset. */
    subset = Calloc1D(cursize + nf, sizeof(int));
    /* initialize the first subset. */
    first_subset(subset + nf, cursize, nf);
    for (i = 0; i < nf; i++)
      subset[i] = i;

    /* iterate over subsets. */
    do {

      /* prepare the current subset. */
      ddata_subset_columns(&dtz, &sub, subset, cursize + nf);
      /* construct the parents' configurations. */
      c_fast_config(sub.col, sub.m.nobs, cursize + nf, sub.nlvl, zptr, &llz, 1);

      c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, nperms, &statistic,
        &pvalue, threshold, type, &df);
      update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ", dtx.m.names[0],
          (pvalue > a) ? "independent from" : "dependent on", dty.m.names[0]);
        for (i = 0; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names, sub.m.ncols));

        Free1D(subset);
        Free1D(zptr);
        FreeDDT(sub);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf, cursize, dtz.m.ncols - nf, nf));

    Free1D(subset);

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_DPERM*/

/* continuous permutation tests. */
static SEXP ast_gperm(gdata dt, int nf, int minsize, int maxsize, double a,
    test_e type, int nperms, double threshold, bool complete, bool debugging) {

int i = 0, j = 0, k = 0, cursize = 0, nc = 0, *subset = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
double **complete_column = NULL;
bool *missing_xy = NULL, *missing_z = NULL;
SEXP retval;
gdata sub = { 0 }, sub_complete = { 0 };

  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));

  if (!complete) {

    missing_xy = Calloc1D(dt.m.nobs, sizeof(bool));
    missing_z = Calloc1D(dt.m.nobs, sizeof(bool));

    gdata_incomplete_cases_range(&dt, missing_xy, 0, 1);

  }/*THEN*/

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset indexes array. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    if (!complete)
      sub_complete = new_gdata(dt.m.nobs, cursize + nf + 2);

    /* iterate over subsets. */
    do {

      /* prepare the current subset. */
      gdata_subset_columns(&dt, &sub, subset, cursize + nf + 2);

      if (!complete) {

        memset(missing_z, '\0', sizeof(bool) * sub.m.nobs);
        gdata_incomplete_cases(&sub, missing_z, 2);
        complete_column = sub_complete.col;

        for (k = 0, nc = 0; k < sub.m.nobs; k++) {

          if (missing_xy[k] || missing_z[k])
            continue;

          for (j = 0; j < sub.m.ncols; j++)
            complete_column[j][nc] = sub.col[j][k];
          nc++;

        }/*FOR*/

      }/*THEN*/
      else {

        complete_column = sub.col;
        nc = sub.m.nobs;

      }/*ELSE*/

      c_gauss_cmcarlo(complete_column, sub.m.ncols, nc, 0, 1, nperms,
        &statistic, &pvalue, threshold, type);
      update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

      /* increment the test counter. */
      test_counter++;

      if (debugging) {

        Rprintf("    > node %s is %s %s given ",
          sub.m.names[0], (pvalue > a) ? "independent from" : "dependent on",
          sub.m.names[1]);
        for (i = 2; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names + 2, sub.m.ncols - 2));

        if (!complete) {

          Free1D(missing_xy);
          Free1D(missing_z);

        }/*THEN*/

        Free1D(subset);
        FreeGDT(sub);
        FreeGDT(sub_complete);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    FreeGDT(sub_complete);
    Free1D(subset);

  }/*FOR*/

  if (!complete) {

    Free1D(missing_xy);
    Free1D(missing_z);

  }/*THEN*/

  FreeGDT(sub);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_GPERM*/

SEXP allsubs_test(SEXP x, SEXP y, SEXP sx, SEXP fixed, SEXP data, SEXP test,
    SEXP B, SEXP alpha, SEXP min, SEXP max, SEXP complete, SEXP debug) {

int minsize = INT(min), maxsize = INT(max);
int i = 0, nf = length(fixed);
double pvalue = 0, min_pvalue = 1, max_pvalue = 0, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
bool debugging = isTRUE(debug);
test_e test_type = test_to_enum(t);
SEXP xx, yy, zz, cc, res = R_NilValue;

  /* call indep_test to deal with tests which have no moving parts, that is when
   * the conditioning set is empty of completely fixed. */
  if (minsize == 0) {

    pvalue = NUM(indep_test(x, y, fixed, data, test, B, alpha, TRUESEXP, complete));
    update_pvalue_range(pvalue, &min_pvalue, &max_pvalue);

    /* increment the test counter. */
    test_counter++;

    if (debugging) {

      Rprintf("    > node %s is %s %s %s",
        CHAR(STRING_ELT(x, 0)),
        (pvalue > a) ? "independent from" : "dependent on",
        CHAR(STRING_ELT(y, 0)),
        (nf > 0) ? "given " : "");
      for (i = 0; i < nf; i++)
        Rprintf("%s ", CHAR(STRING_ELT(fixed, i)));
      Rprintf("(p-value: %g).\n", pvalue);

    }/*THEN*/

    if (pvalue > a) {

      PROTECT(res = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                         a, NULL, 0));
      setAttrib(res, BN_DsepsetSymbol, fixed);

      UNPROTECT(1);
      return res;

    }/*THEN*/
    else {

      /* return even if the variables are not found to be independent, since
       * there are no more tests left to do. */
      if (maxsize == 0)
        return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

    }/*ELSE*/

  }/*THEN*/

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, TRUE));
  PROTECT(yy = c_dataframe_column(data, y, FALSE, TRUE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, TRUE));

  /* extract the missing values indicators. */
  PROTECT(cc = subset_by_name(complete, 3, y, x, sx));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    ddata dtx = ddata_from_SEXP(xx, 0), dty = ddata_from_SEXP(yy, 0);
    ddata dtz = ddata_from_SEXP(zz, 0);
    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dty.m), 0, yy);
    meta_copy_names(&(dtz.m), 0, zz);

    res = ast_discrete(dtx, dty, dtz, nf, minsize, maxsize, test_type, a,
            debugging);

    FreeDDT(dtx);
    FreeDDT(dty);
    FreeDDT(dtz);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    gdata dt = gdata_from_SEXP(zz, 2);
    meta_copy_names(&(dt.m), 2, zz);
    dt.col[0] = REAL(VECTOR_ELT(xx, 0));
    dt.col[1] = REAL(VECTOR_ELT(yy, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));
    dt.m.names[1] = CHAR(STRING_ELT(y, 0));

    if (all_equal(cc, TRUESEXP)) {

      gdata_cache_means(&dt, 0);
      res = ast_gaustests_complete(dt, nf, minsize, maxsize, a, debugging,
              test_type);

    }/*THEN*/
    else {

      res = ast_gaustests_with_missing(dt, nf, minsize, maxsize, a, debugging,
              test_type);

    }/*ELSE*/

    FreeGDT(dt);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    cgdata dtx = cgdata_from_SEXP(xx, 0, 0), dty = cgdata_from_SEXP(yy, 0, 0);
    cgdata dtz = cgdata_from_SEXP(zz, 1, 1);
    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dty.m), 0, yy);
    meta_copy_names(&(dtz.m), 2, zz);

    if (all_equal(cc, TRUESEXP)) {

      res = ast_micg_complete(dtx, dty, dtz, nf, minsize, maxsize, a,
              debugging);

    }/*THEN*/
    else {

      res = ast_micg_with_missing(dtx, dty, dtz, nf, minsize, maxsize, a,
              debugging);

    }/*ELSE*/

    FreeCGDT(dtx);
    FreeCGDT(dty);
    FreeCGDT(dtz);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    ddata dtx = ddata_from_SEXP(xx, 0), dty = ddata_from_SEXP(yy, 0);
    ddata dtz = ddata_from_SEXP(zz, 0);
    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dty.m), 0, yy);
    meta_copy_names(&(dtz.m), 0, zz);

    res = ast_dperm(dtx, dty, dtz, nf, minsize, maxsize, a, test_type,
            INT(B), IS_SMC(test_type) ? a : 1, debugging);

    FreeDDT(dtx);
    FreeDDT(dty);
    FreeDDT(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    gdata dt = gdata_from_SEXP(zz, 2);
    meta_copy_names(&(dt.m), 2, zz);
    dt.col[0] = REAL(VECTOR_ELT(xx, 0));
    dt.col[1] = REAL(VECTOR_ELT(yy, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));
    dt.m.names[1] = CHAR(STRING_ELT(y, 0));
    gdata_cache_means(&dt, 0);

    res = ast_gperm(dt, nf, minsize, maxsize, a, test_type,
            INT(B), IS_SMC(test_type) ? a : 1, all_equal(cc, TRUESEXP),
            debugging);

    FreeGDT(dt);

  }/*THEN*/

  /* the p-values of the tests performed under minsize == 0 should be considered
   * when computing the maximum and minimum p-values. */
  if (minsize == 0) {

    REAL(res)[1] = REAL(res)[1] < min_pvalue ? REAL(res)[1] : min_pvalue;
    REAL(res)[2] = REAL(res)[2] > max_pvalue ? REAL(res)[2] : max_pvalue;

  }/*THEN*/

  UNPROTECT(4);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  return res;

}/*ALLSUBS_TEST*/

