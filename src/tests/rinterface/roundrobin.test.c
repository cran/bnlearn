#include "../../include/rcore.h"
#include "../../core/data.table.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../patterns.h"
#include "../tests.h"

SEXP roundrobin_test(SEXP x, SEXP z, SEXP fixed, SEXP data, SEXP test,
    SEXP alpha, SEXP extra_args, SEXP complete, SEXP debug) {

double *pvalue = NULL, a = NUM(alpha);
bool debugging = isTRUE(debug);
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_to_enum(t);
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
    tabular dtx = tabular_from_SEXP(xx, 0, 0), dtz = tabular_from_SEXP(zz, 0, 0);

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 0, zz);
    meta_init_flags(&(dtz.m), 0, R_NilValue, which_fixed);

    rrd_discrete(dtx, dtz, test_type, pvalue, a, debugging);

    FreeTAB(dtx);
    FreeTAB(dtz);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    tabular dt = tabular_from_SEXP(zz, 0, 1);

    meta_copy_names(&(dt.m), 1, zz);
    meta_init_flags(&(dt.m), 1, R_NilValue, which_fixed);
    dt.ccol[0] = REAL(VECTOR_ELT(xx, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));

    if (all_equal(cc, TRUESEXP)) {

      tabular_cache_means(&dt, 0);
      rrd_gaustests_complete(dt, test_type, pvalue, a, debugging);

    }/*THEN*/
    else {

      rrd_gaustests_with_missing(dt, test_type, pvalue, a, debugging);

    }/*ELSE*/

    FreeTAB(dt);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian test. */
    tabular dtx = tabular_from_SEXP(xx, 0, 0), dtz = tabular_from_SEXP(zz, 1, 1);

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 2, zz);
    meta_init_flags(&(dtz.m), 2, R_NilValue, which_fixed);

    if (all_equal(cc, TRUESEXP)) {

      rrd_micg_complete(dtx, dtz, test_type, pvalue, a, debugging);

    }/*THEN*/
    else {

      rrd_micg_with_missing(dtx, dtz, test_type, pvalue, a, debugging);

    }/*ELSE*/

    FreeTAB(dtx);
    FreeTAB(dtz);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    tabular dtx = tabular_from_SEXP(xx, 0, 0), dtz = tabular_from_SEXP(zz, 0, 0);
    int B = INT(getListElement(extra_args, "B"));

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 0, zz);
    meta_init_flags(&(dtz.m), 0, R_NilValue, which_fixed);

    rrd_dperm(dtx, dtz, test_type, pvalue, a, B, IS_SMC(test_type) ? a : 1,
      debugging);

    FreeTAB(dtx);
    FreeTAB(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    tabular dt = tabular_from_SEXP(zz, 0, 1);
    int B = INT(getListElement(extra_args, "B"));

    meta_copy_names(&(dt.m), 1, zz);
    meta_init_flags(&(dt.m), 1, R_NilValue, which_fixed);
    dt.ccol[0] = REAL(VECTOR_ELT(xx, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));

    rrd_gperm(dt, test_type, pvalue, a, B, IS_SMC(test_type) ? a : 1,
      all_equal(cc, TRUESEXP), debugging);

    FreeTAB(dt);

  }/*THEN*/
  else if (test_type == CUSTOM_T) {

    /* user-provided test function. */
    SEXP custom_fn = getListElement(extra_args, "fun");
    SEXP custom_args = getListElement(extra_args, "args");

    rrd_custom(x, z, fixed, data, custom_fn, custom_args, pvalue, a, debugging);

  }/*THEN*/

  UNPROTECT(5);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  /* increment the test counter. */
  test_counter += length(zz) - length(fixed);

  return result;

}/*ROUNDROBIN_TEST*/
