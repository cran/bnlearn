#include "include/rcore.h"

/* return the unique elements from an input vector.*/
SEXP unique(SEXP array) {

int *d = NULL, i = 0, k = 0, dup_counter = 0, n = length(array);
int *res = NULL, *a = NULL;
SEXP dup, result = R_NilValue;

  PROTECT(dup = duplicated(array, FALSE));
  d = LOGICAL(dup);

  switch(TYPEOF(array)) {

    case INTSXP:

      a = INTEGER(array);

      for (i = 0; i < n; i++)
        if ((d[i] == 0) && (a[i] != NA_INTEGER))
          dup_counter++;

      PROTECT(result = allocVector(INTSXP, dup_counter));
      res = INTEGER(result);

      for (i = 0; i < n; i++)
        if ((d[i] == 0) && (a[i] != NA_INTEGER))
          res[k++] = a[i];

      break;

    case STRSXP:

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          dup_counter++;

      PROTECT(result = allocVector(STRSXP, dup_counter));

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          SET_STRING_ELT(result, k++, STRING_ELT(array, i));

      break;

    default:

      error("this SEXP type is not handled in unique().");

  }/*SWITCH*/

  UNPROTECT(2);

  return result;

}/*UNIQUE*/

/* determine which elements are dupes. */
SEXP dupe(SEXP array) {

int i = 0, n = length(array);
int *res = NULL, *tmp = NULL;
SEXP result, temp;

  PROTECT(result = duplicated(array, FALSE));
  PROTECT(temp = duplicated(array, TRUE));
  res = LOGICAL(result);
  tmp = LOGICAL(temp);

  for (i = 0; i < n; i++)
    res[i] = res[i] || tmp[i];

  UNPROTECT(2);

  return result;

}/*DUPE*/

