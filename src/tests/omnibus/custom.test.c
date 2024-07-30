#include "../../include/rcore.h"
#include "../tests.h"

double custom_test_function(SEXP x, SEXP y, SEXP z, SEXP data, SEXP custom_fn,
    SEXP custom_args, double *pvalue) {

double statistic = 0, temp_pvalue = 0;
SEXP call, args_iterator, result;

  /* allocate and populate the pairlist to be valuated. */
  PROTECT(args_iterator = call = allocLang(6));
  /* first slot, the function name. */
  SETCAR(args_iterator, custom_fn);
  args_iterator = CDR(args_iterator);
  /* second slot, the label of the first node. */
  SETCAR(args_iterator, x);
  args_iterator = CDR(args_iterator);
  /* third slot, the label of the second node. */
  SETCAR(args_iterator, y);
  args_iterator = CDR(args_iterator);
  /* fourth slot, the labels of the conditioning set (NULL if empty). */
  SETCAR(args_iterator, z);
  args_iterator = CDR(args_iterator);
  /* fifth slot, the data. */
  SETCAR(args_iterator, data);
  args_iterator = CDR(args_iterator);
  /* sixth slot, the optional arguments passed as a list. */
  SETCAR(args_iterator, custom_args);
  /* evaluate the custom score function. */
  PROTECT(result = eval(call, R_GlobalEnv));

  /* the return value must be a scalar, real number. */
  if ((TYPEOF(result) != REALSXP) || (length(result) != 2))
    error("the test for nodes %s and %s must return two scalar, real values.",
      CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)));

  /* the second element of the return value is a p-value. */
  temp_pvalue = REAL(result)[1];
  if (ISNAN(temp_pvalue))
    error("the test for nodes %s and %s has a NA p-value.",
        CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)));
  if ((temp_pvalue < 0) || (temp_pvalue > 1))
    error("the test for nodes %s and %s has a p-value not in [0, 1].",
        CHAR(STRING_ELT(x, 0)), CHAR(STRING_ELT(y, 0)));

  statistic = REAL(result)[0];
  *pvalue = temp_pvalue;

  UNPROTECT(2);

  return statistic;

}/*CUSTOM_TEST_FUNCTION*/


