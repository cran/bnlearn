#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/sets.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/common.h"
#include "../../include/globals.h"
#include "../../core/covariance.matrix.h"
#include "../../core/correlation.h"
#include "../../core/data.table.h"
#include "../../math/linear.algebra.h"
#include "../tests.h"
#include "../patterns.h"

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

