#include "include/rcore.h"

/* setdiff() for two vectors of character strings. */
SEXP string_setdiff(SEXP large, SEXP small) {

int i = 0, k = 0, *t = NULL, nl = length(large), ns = length(small);
SEXP try, diff;

  /* match the elements of the smaller set against the larger set. */
  PROTECT(try = match(small, large, 0));
  t = INTEGER(try);
  /* allocate the return value. */
  PROTECT(diff = allocVector(STRSXP, nl - ns));
  /* copy the elements that are not shared between te two. */
  for (i = 0, k = 0; i < nl; i++)
    if (t[i] == 0)
      SET_STRING_ELT(diff, k++, STRING_ELT(large, i));

  UNPROTECT(2);

  return diff;

}/*STRING_SETDIFF*/

/* delete a string for a STRSXP array. */
SEXP string_delete(SEXP array, SEXP string, int *idx) {

int i = 0, k = 0, *t = NULL, n = length(array);
SEXP try, new_array;

  PROTECT(try = match(array, string, 0));
  t = INTEGER(try);

  /* optional, save the index of the string. */
  if (idx)
    *idx = *t;

  /*allocate the new arra. */
  PROTECT(new_array = allocVector(STRSXP, n - 1));

  for (i = 0, k = 0; i < n; i++)
    if (i != *t - 1)
      SET_STRING_ELT(new_array, k++, STRING_ELT(array, i));

  UNPROTECT(2);

  return new_array;

}/*STRING_DELETE*/

/* variadic version of mkString(). */
SEXP mkStringVec(int n, ...) {

va_list strings;
int i = 0;
SEXP vec;

  PROTECT(vec = allocVector(STRSXP, n));
  va_start(strings, n);
  for (i = 0; i < n; i++)
    SET_STRING_ELT(vec, i, mkChar(va_arg(strings, char *)));
  va_end(strings);
  UNPROTECT(1);

  return vec;

}/*MKSTRINGVEC*/
