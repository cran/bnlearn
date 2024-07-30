#include "../../include/rcore.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../minimal/common.h"
#include "../../include/globals.h"
#include "../../core/data.table.h"
#include "../tests.h"
#include "../patterns.h"

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
    int B = INT(getListElement(extra_args, "B"));

    meta_copy_names(&(dtx.m), 0, xx);
    meta_copy_names(&(dtz.m), 0, zz);
    meta_init_flags(&(dtz.m), 0, R_NilValue, which_fixed);

    rrd_dperm(dtx, dtz, test_type, pvalue, a, B, IS_SMC(test_type) ? a : 1,
      debugging);

    FreeDDT(dtx);
    FreeDDT(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    gdata dt = gdata_from_SEXP(zz, 1);
    int B = INT(getListElement(extra_args, "B"));

    meta_copy_names(&(dt.m), 1, zz);
    meta_init_flags(&(dt.m), 1, R_NilValue, which_fixed);
    dt.col[0] = REAL(VECTOR_ELT(xx, 0));
    dt.m.names[0] = CHAR(STRING_ELT(x, 0));

    rrd_gperm(dt, test_type, pvalue, a, B, IS_SMC(test_type) ? a : 1,
      all_equal(cc, TRUESEXP), debugging);

    FreeGDT(dt);

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
