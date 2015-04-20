#include "include/rcore.h"
#include "include/tests.h"

/* independence tests, frontend to be used in R code. */
SEXP indep_test(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP learning) {

  /* if either node to test is not provided, return a zero-length numeric
   * vector. */
  if (length(x) == 0 || length(y) == 0)
    return allocVector(REALSXP, 0);

  /* filter for NULL and empty strings to make it easy to interface with R. */
  if (length(sx) == 0 || sx == R_NilValue)
    return utest(x, y, data, test, B, alpha, learning);
  else
    return ctest(x, y, sx, data, test, B, alpha, learning);

}/*INDEP_TEST*/
