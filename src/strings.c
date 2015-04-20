#include "include/rcore.h"

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
