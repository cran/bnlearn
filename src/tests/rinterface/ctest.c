#include "../../include/rcore.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/common.h"
#include "../../include/globals.h"
#include "../../core/data.table.h"
#include "../tests.h"
#include "../patterns.h"

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
    ddata dtx = ddata_from_SEXP(xx2, 0), dty = ddata_from_SEXP(yy2, 0);
    ddata dtz = ddata_from_SEXP(zz, 0);

    statistic = ct_discrete(dtx, dty, dtz, pvalue, &df, test_type);

    FreeDDT(dtx);
    FreeDDT(dty);
    FreeDDT(dtz);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    gdata dtx = gdata_from_SEXP(xx2, 0), dt = gdata_from_SEXP(zz, 2);
    dt.col[1] = REAL(VECTOR_ELT(yy2, 0));

    if (all_equal(cc, TRUESEXP)) {

      gdata_cache_means(&dt, 1);
      statistic = ct_gaustests_complete(dtx, dt, pvalue, &df, test_type);

    }/*THEN*/
    else {

      statistic = ct_gaustests_with_missing(dtx, dt, pvalue, &df, test_type);

    }/*ELSE*/

    FreeGDT(dtx);
    FreeGDT(dt);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    cgdata dtx = cgdata_from_SEXP(xx2, 0, 0), dty = cgdata_from_SEXP(yy2, 0, 0);
    cgdata dtz = cgdata_from_SEXP(zz, 1, 1);

    if (all_equal(cc, TRUESEXP))
      statistic = ct_micg_complete(dtx, dty, dtz, pvalue, &df);
    else
      statistic = ct_micg_with_missing(dtx, dty, dtz, pvalue, &df);

    FreeCGDT(dtx);
    FreeCGDT(dty);
    FreeCGDT(dtz);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    ddata dtx = ddata_from_SEXP(xx2, 0), dty = ddata_from_SEXP(yy2, 0);
    ddata dtz = ddata_from_SEXP(zz, 0);

    int B = INT(getListElement(extra_args, "B"));

    statistic = ct_dperm(dtx, dty, dtz, pvalue, &df, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1);

    FreeDDT(dtx);
    FreeDDT(dty);
    FreeDDT(dtz);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    gdata dtx = gdata_from_SEXP(xx2, 0), dt = gdata_from_SEXP(zz, 2);
    dt.col[1] = REAL(VECTOR_ELT(yy2, 0));

    int B = INT(getListElement(extra_args, "B"));

    statistic = ct_gperm(dtx, dt, pvalue, &df, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1, all_equal(cc, TRUESEXP));

    FreeGDT(dtx);
    FreeGDT(dt);

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

