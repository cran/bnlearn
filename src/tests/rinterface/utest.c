#include "../../include/rcore.h"
#include "../../include/globals.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/common.h"
#include "../tests.h"
#include "../patterns.h"

/* unconditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP alpha, SEXP extra_args,
    SEXP learning, SEXP complete) {

int ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
test_e test_type = test_to_enum(t);
SEXP xx, yy, cc, result;

  /* allocate the return value, which has the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  nobs = length(yy);

  /* extract the missing values indicators. */
  PROTECT(cc = subset_by_name(complete, 2, y, x));

  if (IS_DISCRETE_ASYMPTOTIC_TEST(test_type)) {

    /* parametric tests for discrete variables. */
    statistic = ut_discrete(xx, yy, nobs, ntests, pvalue, &df, test_type);

  }/*THEN*/
  else if ((test_type == COR) || (test_type == ZF) || (test_type == MI_G) ||
           (test_type == MI_G_SH)) {

    /* parametric tests for Gaussian variables. */
    if (all_equal(cc, TRUESEXP))
      statistic = ut_gaustests_complete(xx, yy, nobs, ntests, pvalue,
                    &df, test_type);
    else
      statistic = ut_gaustests_with_missing(xx, yy, nobs, ntests, pvalue,
                    &df, test_type);

  }/*THEN*/
  else if (test_type == MI_CG) {

    /* conditional linear Gaussian mutual information test. */
    if (all_equal(cc, TRUESEXP))
      statistic = ut_micg_complete(xx, yy, nobs, ntests, pvalue, &df);
    else
      statistic = ut_micg_with_missing(xx, yy, nobs, ntests, pvalue, &df);

  }/*THEN*/
  else if (IS_DISCRETE_PERMUTATION_TEST(test_type)) {

    /* discrete permutation tests. */
    int B = INT(getListElement(extra_args, "B"));

    statistic = ut_dperm(xx, yy, nobs, ntests, pvalue, &df, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1);

  }/*THEN*/
  else if (IS_CONTINUOUS_PERMUTATION_TEST(test_type)) {

    /* continuous permutation tests. */
    int B = INT(getListElement(extra_args, "B"));

    statistic = ut_gperm(xx, yy, nobs, ntests, pvalue, test_type, B,
                  IS_SMC(test_type) ? NUM(alpha) : 1, all_equal(cc, TRUESEXP));

  }/*THEN*/
  else if (test_type == CUSTOM_T) {

    /* user-provided test function. */
    SEXP custom_fn = getListElement(extra_args, "fun");
    SEXP custom_args = getListElement(extra_args, "args");

    statistic = ut_custom(x, y, data, custom_fn, custom_args, pvalue);

  }/*THEN*/

  UNPROTECT(4);

  /* catch-all for unknown tests (after deallocating memory.) */
  if (test_type == ENOTEST)
    error("unknown test statistic '%s'.", t);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, extra_args);

}/*UTEST*/

