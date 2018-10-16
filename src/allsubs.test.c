#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/data.table.h"

#define PVALUE(testfun) \
        pvalue = testfun; \
        min_pvalue = pvalue < min_pvalue ? pvalue : min_pvalue; \
        max_pvalue = pvalue > max_pvalue ? pvalue : max_pvalue; \
        /* increment the test counter. */ \
        test_counter++;

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
static SEXP ast_discrete(SEXP xx, SEXP yy, SEXP zz, SEXP x, SEXP y, int nf,
    int minsize, int maxsize, test_e test, double a, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
double statistic = 0, pvalue = 0, df = 0, min_pvalue = 1, max_pvalue = 0;
SEXP retval;
ddata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = ddata_from_SEXP(zz, 0);
  meta_copy_names(&(dt.m), 0, zz);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dt.m.nobs, sizeof(int));

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
      ddata_subset_columns(&dt, &sub, subset, cursize + nf);
      /* construct the parents' configurations. */
      c_fast_config(sub.col, sub.m.nobs, cursize + nf, sub.nlvl, zptr, &llz, 1);

      if (test == MI || test == MI_ADF || test == X2 || test == X2_ADF) {

        /* mutual information and Pearson's X^2 asymptotic tests. */
        statistic = c_cchisqtest(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, &df,
                      test, (test == MI) || (test == MI_ADF));
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_SH) {

        /* shrinkage mutual information test. */
        statistic = c_shcmi(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, &df, TRUE);
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == JT) {

        /* Jonckheere-Terpstra test. */
        statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs);
        PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));

      }/*THEN*/

      if (debuglevel > 0) {

        Rprintf("    > node %s is %s %s given ",
          CHAR(STRING_ELT(x, 0)),
          (pvalue > a) ? "independent from" : "dependent on",
          CHAR(STRING_ELT(y, 0)));
        for (i = 0; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names, sub.m.ncols));

        Free1D(subset);
        Free1D(zptr);
        FreeDDT(dt, FALSE);
        FreeDDT(sub, FALSE);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf, cursize, dt.m.ncols - nf, nf));

    Free1D(subset);

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(dt, FALSE);
  FreeDDT(sub, FALSE);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_DISCRETE*/

/* parametric tests for Gaussian variables. */
static SEXP ast_gaustests(SEXP xx, SEXP yy, SEXP zz, SEXP x, SEXP y, int nf,
    int minsize, int maxsize, double a, int debuglevel, test_e test) {

int i = 0, cursize = 0, *subset = NULL, df = 0;
double statistic = 0, lambda = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
SEXP retval;
gdata dt = { 0 }, sub = { 0 };
covariance cov = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 2);
  meta_copy_names(&(dt.m), 2, zz);
  dt.col[0] = REAL(xx);
  dt.col[1] = REAL(yy);
  dt.m.names[0] = CHAR(STRING_ELT(x, 0));
  dt.m.names[1] = CHAR(STRING_ELT(y, 0));
  /* allocate and compute mean values. */
  gdata_cache_means(&dt, 0);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* compute the degrees of freedom for correlation and mutual information. */
    if (test == COR)
      df = dt.m.nobs - (cursize + nf + 2);
    else if ((test == MI_G) || (test == MI_G_SH))
      df = 1;

    if (((test == COR) && (df < 1)) ||
        ((test == ZF) && (dt.m.nobs - (cursize + nf + 2) < 2))) {

      /* if there are not enough degrees of freedom, return independence. */
      warning("trying to do a conditional independence test with zero degrees of freedom.");

      PVALUE(1);

      if (debuglevel > 0) {

        Rprintf("    > node %s is independent from %s given any conditioning set of size %d ",
          CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)), cursize + nf);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      /* the conditioning set is the first possible conditioning set, which
       * comprises the first columns of the data table. */
      PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                         a, dt.m.names + 2, cursize + nf));

      FreeGDT(dt, FALSE);
      FreeGDT(sub, FALSE);

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
        PVALUE(2 * pt(fabs(statistic), df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_G) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == MI_G_SH) {

        lambda = covmat_lambda(sub.col, sub.mean, cov, sub.m.nobs, NULL,
                   sub.m.nobs);
        covmat_shrink(cov, lambda);
        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = 2 * sub.m.nobs * cor_mi_trans(statistic);
        PVALUE(pchisq(statistic, df, FALSE, FALSE));

      }/*THEN*/
      else if (test == ZF) {

        statistic = c_fast_pcor(cov, 0, 1, NULL, TRUE);
        statistic = cor_zf_trans(statistic, (double)sub.m.nobs - sub.m.ncols);
        PVALUE(2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE));

      }/*THEN*/

      if (debuglevel > 0) {

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
        FreeGDT(dt, FALSE);
        FreeGDT(sub, FALSE);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    FreeCOV(cov);
    Free1D(subset);

  }/*FOR*/

  FreeGDT(dt, FALSE);
  FreeGDT(sub, FALSE);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_GAUSTESTS*/

/* conditional linear Gaussian test. */
static SEXP ast_micg(SEXP xx, SEXP yy, SEXP zz, SEXP x, SEXP y, int nf,
    int minsize, int maxsize, double a, int debuglevel) {

int cursize = 0, xtype = TYPEOF(xx), ytype = TYPEOF(yy);
int i = 0, *subset = NULL;
int *zptr = NULL, llx = 0, lly = 0, llz = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
void *xptr = 0, *yptr = 0;
SEXP retval;
cgdata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = cgdata_from_SEXP(zz, 1, 1);
  meta_copy_names(&(dt.m), 2, zz);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_cgdata(dt.m.nobs, dt.ndcols, dt.ngcols);

  /* if both variables are continuous and all conditioning variables are
   * continuous, the test reverts back to a Gaussian mutual information test. */
  if ((xtype == INTSXP) && (ytype == INTSXP)) {

    xptr = INTEGER(xx);
    llx = NLEVELS(xx);
    yptr = INTEGER(yy);
    lly = NLEVELS(yy);

  }/*THEN*/
  if ((xtype == REALSXP) && (ytype == REALSXP)) {

    xptr = REAL(xx);
    yptr = REAL(yy);

  }/*THEN*/
  else if ((xtype == REALSXP) && (ytype == INTSXP)) {

    xptr = REAL(xx);
    yptr = INTEGER(yy);
    lly = NLEVELS(yy);

  }/*THEN*/
  else if ((xtype == INTSXP) && (ytype == REALSXP)) {

    yptr = INTEGER(xx);
    lly = NLEVELS(xx);
    xptr = REAL(yy);

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
      cgdata_subset_columns(&dt, &sub, subset, cursize + nf + 2);

      /* if there are discrete conditioning variables, compute their
       * configurations. */
      if (sub.ndcols - 1 > 0) {

        zptr = Calloc1D(sub.m.nobs, sizeof(int));
        c_fast_config(sub.dcol + 1, sub.m.nobs, sub.ndcols - 1, sub.nlvl + 1,
          zptr, &llz, 1);

      }/*THEN*/

      if ((ytype == INTSXP) && (xtype == INTSXP)) {

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
      else if ((ytype == REALSXP) && (xtype == REALSXP)) {

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

      PVALUE(pchisq(statistic, df, FALSE, FALSE));

      if (debuglevel > 0) {

        Rprintf("    > node %s is %s %s given ",
          CHAR(STRING_ELT(x, 0)),
          (pvalue > a) ? "independent from" : "dependent on",
          CHAR(STRING_ELT(y, 0)));
        for (i = 2; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names + 2, sub.m.ncols - 2));

        Free1D(subset);
        FreeCGDT(dt, FALSE);
        FreeCGDT(sub, FALSE);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    Free1D(subset);

  }/*FOR*/

  FreeCGDT(dt, FALSE);
  FreeCGDT(sub, FALSE);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_MICG*/

/* discrete permutation tests. */
static SEXP ast_dperm(SEXP xx, SEXP yy, SEXP zz, SEXP x, SEXP y, int nf,
    int minsize, int maxsize, double a, test_e type, int nperms,
    double threshold, int debuglevel) {

int *xptr = INTEGER(xx), *yptr = INTEGER(yy), *zptr = NULL, *subset = NULL;
int i = 0, cursize = 0, llx = NLEVELS(xx), lly = NLEVELS(yy), llz = 0;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0, df = 0;
SEXP retval;
ddata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = ddata_from_SEXP(zz, 0);
  meta_copy_names(&(dt.m), 0, zz);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* allocate the parents' configurations. */
  zptr = Calloc1D(dt.m.nobs, sizeof(int));

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
      ddata_subset_columns(&dt, &sub, subset, cursize + nf);
      /* construct the parents' configurations. */
      c_fast_config(sub.col, sub.m.nobs, cursize + nf, sub.nlvl, zptr, &llz, 1);

      c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, sub.m.nobs, nperms, &statistic,
        &pvalue, threshold, type, &df);
      PVALUE(pvalue);

      if (debuglevel > 0) {

        Rprintf("    > node %s is %s %s given ",
          CHAR(STRING_ELT(x, 0)),
          (pvalue > a) ? "independent from" : "dependent on",
          CHAR(STRING_ELT(y, 0)));
        for (i = 0; i < sub.m.ncols; i++)
          Rprintf("%s ", sub.m.names[i]);
        Rprintf("(p-value: %g).\n", pvalue);

      }/*THEN*/

      if (pvalue > a) {

        PROTECT(retval = ast_prepare_retval(pvalue, min_pvalue, max_pvalue,
                           a, sub.m.names, sub.m.ncols));

        Free1D(subset);
        Free1D(zptr);
        FreeDDT(dt, FALSE);
        FreeDDT(sub, FALSE);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf, cursize, dt.m.ncols - nf, nf));

    Free1D(subset);

  }/*FOR*/

  Free1D(zptr);
  FreeDDT(dt, FALSE);
  FreeDDT(sub, FALSE);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_DPERM*/

/* continuous permutation tests. */
static SEXP ast_gperm(SEXP xx, SEXP yy, SEXP zz, SEXP x, SEXP y, int nf,
    int minsize, int maxsize, double a, test_e type, int nperms,
    double threshold, int debuglevel) {

int i = 0, cursize = 0, *subset = NULL;
double statistic = 0, pvalue = 0, min_pvalue = 1, max_pvalue = 0;
SEXP retval;
gdata dt = { 0 }, sub = { 0 };

  /* allocate and initialize a data table for the variables. */
  dt = gdata_from_SEXP(zz, 2);
  meta_copy_names(&(dt.m), 2, zz);
  dt.col[0] = REAL(xx);
  dt.col[1] = REAL(yy);
  dt.m.names[0] = CHAR(STRING_ELT(x, 0));
  dt.m.names[1] = CHAR(STRING_ELT(y, 0));
  /* allocate and compute mean values. */
  gdata_cache_means(&dt, 0);
  /* allocate a second data table to hold the conditioning variables. */
  sub = empty_gdata(dt.m.nobs, dt.m.ncols);
  sub.mean = Calloc1D(dt.m.ncols, sizeof(double));

  for (cursize = fmax(1, minsize); cursize <= maxsize; cursize++) {

    /* allocate and initialize the subset indexes array. */
    subset = Calloc1D(cursize + nf + 2, sizeof(int));
    /* initialize the first subset. */
    first_subset(subset + nf + 2, cursize, nf + 2);
    for (i = 0; i < nf + 2; i++)
      subset[i] = i;

    /* iterate over subsets. */
    do {

      /* prepare the current subset. */
      gdata_subset_columns(&dt, &sub, subset, cursize + nf + 2);

      c_gauss_cmcarlo(sub.col, sub.m.ncols, sub.m.nobs, 0, 1, nperms,
        &statistic, &pvalue, threshold, type);
      PVALUE(pvalue);

      if (debuglevel > 0) {

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
        FreeGDT(dt, FALSE);
        FreeGDT(sub, FALSE);

        UNPROTECT(1);
        return retval;

      }/*THEN*/

    } while (next_subset(subset + nf + 2, cursize, dt.m.ncols - nf - 2, nf + 2));

    Free1D(subset);

  }/*FOR*/

  FreeGDT(dt, FALSE);
  FreeGDT(sub, FALSE);

  return ast_prepare_retval(pvalue, min_pvalue, max_pvalue, a, NULL, 0);

}/*AST_GPERM*/

SEXP allsubs_test(SEXP x, SEXP y, SEXP sx, SEXP fixed, SEXP data, SEXP test,
    SEXP B, SEXP alpha, SEXP min, SEXP max, SEXP complete, SEXP debug) {

int minsize = INT(min), maxsize = INT(max), debuglevel = isTRUE(debug);
int i = 0, nf = length(fixed);
double pvalue = 0, min_pvalue = 1, max_pvalue = 0, a = NUM(alpha);
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_label(t);
SEXP xx, yy, zz, res = R_NilValue;

  /* call indep_test to deal with zero-length conditioning subsets. */
  if (minsize == 0) {

    PVALUE(NUM(indep_test(x, y, fixed, data, test, B, alpha, TRUESEXP, complete)));

    if (debuglevel > 0) {

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
  PROTECT(xx = c_dataframe_column(data, x, TRUE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, TRUE));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    res = ast_discrete(xx, yy, zz, x, y, nf, minsize, maxsize, test_type, a,
            debuglevel);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    res = ast_gaustests(xx, yy, zz, x, y, nf, minsize, maxsize, a,
            debuglevel, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    res = ast_micg(xx, yy, zz, x, y, nf, minsize, maxsize, a, debuglevel);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    res = ast_dperm(xx, yy, zz, x, y, nf, minsize, maxsize, a, test_type,
            INT(B), IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    res = ast_gperm(xx, yy, zz, x, y, nf, minsize, maxsize, a, test_type,
            INT(B), IS_SMC(test_type) ? a : 1, debuglevel);

  }/*THEN*/

  UNPROTECT(3);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  return res;

}/*ALLSUBS_TEST*/

