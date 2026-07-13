#include "../../include/rcore.h"
#include "../../core/data.table.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/data.frame.h"
#include "../patterns.h"
#include "../tests.h"

/* conditional independence tests. */
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP alpha,
    SEXP extra_args, SEXP learning, SEXP complete) {

int ntests = length(x);
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_to_enum(t);
SEXP xx2, yy2, zz, cc, result;

  /* allocate the return value, which the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx2 = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy2 = c_dataframe_column(data, y, FALSE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, FALSE));

  /* extract the missing values indicators. */
  PROTECT(cc = subset_by_name(complete, 3, y, x, sx));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    tabular dtx = tabular_from_SEXP(xx2, 0, 0), dty = tabular_from_SEXP(yy2, 0, 0);
    tabular dtz = tabular_from_SEXP(zz, 0, 0);

    statistic = ct_discrete(dtx, dty, dtz, pvalue, &df, test_type);

    FreeTAB(dtx);
    FreeTAB(dty);
    FreeTAB(dtz);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    tabular dtx = tabular_from_SEXP(xx2, 0, 0), dt = tabular_from_SEXP(zz, 0, 2);
    dt.ccol[1] = REAL(VECTOR_ELT(yy2, 0));

    if (all_equal(cc, TRUESEXP)) {

      tabular_cache_means(&dt, 1);
      statistic = ct_gaustests_complete(dtx, dt, pvalue, &df, test_type);

    }/*THEN*/
    else {

      statistic = ct_gaustests_with_missing(dtx, dt, pvalue, &df, test_type);

    }/*ELSE*/

    FreeTAB(dtx);
    FreeTAB(dt);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    tabular dtx = tabular_from_SEXP(xx2, 0, 0), dty = tabular_from_SEXP(yy2, 0, 0);
    tabular dtz = tabular_from_SEXP(zz, 1, 1);

    if (all_equal(cc, TRUESEXP))
      statistic = ct_micg_complete(dtx, dty, dtz, pvalue, &df);
    else
      statistic = ct_micg_with_missing(dtx, dty, dtz, pvalue, &df);

    FreeTAB(dtx);
    FreeTAB(dty);
    FreeTAB(dtz);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    tabular dtx = tabular_from_SEXP(xx2, 0, 0), dty = tabular_from_SEXP(yy2, 0, 0);
    tabular dtz = tabular_from_SEXP(zz, 0, 0);

    int B = INT(getListElement(extra_args, "B"));

    statistic = ct_dperm(dtx, dty, dtz, pvalue, &df, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1);

    FreeTAB(dtx);
    FreeTAB(dty);
    FreeTAB(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    tabular dtx = tabular_from_SEXP(xx2, 0, 0), dt = tabular_from_SEXP(zz, 0, 2);
    dt.ccol[1] = REAL(VECTOR_ELT(yy2, 0));

    int B = INT(getListElement(extra_args, "B"));

    statistic = ct_gperm(dtx, dt, pvalue, &df, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1, all_equal(cc, TRUESEXP));

    FreeTAB(dtx);
    FreeTAB(dt);

  }/*THEN*/
  else if (test_type == CUSTOM_T) {

    /* user-provided test function. */
    SEXP custom_fn = getListElement(extra_args, "fun");
    SEXP custom_args = getListElement(extra_args, "args");

    statistic = ct_custom(x, y, sx, data, custom_fn, custom_args, pvalue);

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
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, extra_args);

}/*CTEST*/

