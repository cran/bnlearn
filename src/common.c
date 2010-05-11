#include <R.h>
#include <Rinternals.h>

/* get the list element named str, or return NULL. */
SEXP getListElement(SEXP list, char *str) {

SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
int i = 0;

  for (i = 0; i < length(list); i++) {

    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {

      elmt = VECTOR_ELT(list, i);
      break;

    }/*THEN*/

  }/*FOR*/

return elmt;

}/*GETLISTELEMENT*/

/* sampling without replacement, internal copy of the
 * SampleNoReplace function in src/main/random.c.
 * Copyright (C) 2003--2008  The R Foundation,
 * licensed under "GPLv2 or later" licence. */
void SampleNoReplace(int k, int n, int *y, int *x) {

int i = 0, j = 0;

  for (i = 0; i < n; i++)
    x[i] = i;
  for (i = 0; i < k; i++) {

    j = n * unif_rand();
    y[i] = x[j] + 1;
    x[j] = x[--n];

  }/*FOR*/

}/*SAMPLENOREPLACE*/

/* return the unique elements from an input vector.*/
SEXP unique(SEXP array) {

int *d = NULL, i = 0, k = 0, dup_counter = 0, n = LENGTH(array);
SEXP dup, result = R_NilValue;

  PROTECT(dup = duplicated(array, FALSE));
  d = LOGICAL(dup);

  for (i = 0; i < n; i++)
    if (d[i] == 0)
      dup_counter++;

  switch(TYPEOF(array)) {

    case INTSXP:

      PROTECT(result = allocVector(INTSXP, dup_counter));
      int *res = INTEGER(result), *a = INTEGER(array);

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          res[k++] = a[i];

      break;

    case STRSXP:

      PROTECT(result = allocVector(STRSXP, dup_counter));

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          SET_STRING_ELT(result, k++, STRING_ELT(array, i));

      break;

    default:

      error("this SEXP type is no handled in unique().");

  }/*SWITCH*/

  UNPROTECT(2);

  return result;

}/*UNIQUE*/

/* determine which elements are dupes. */
SEXP dupe(SEXP array) {

int i = 0, n = LENGTH(array);
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

/* transform an integer vector into a factor. */
SEXP int2fac(SEXP vector) {

SEXP result, class, levels, lvls;

  /* find out which levels are needed. */
  PROTECT(levels = unique(vector));

  /* match the elements of the vector against the levels. */
  PROTECT(result = match(levels, vector, 0));

  /* set the levels of the factor. */
  PROTECT(lvls = coerceVector(levels, STRSXP));
  setAttrib(result, R_LevelsSymbol, lvls);
 
  /* set the class of the return value. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(result, R_ClassSymbol, class);

  UNPROTECT(4);

  return result;

}/*INT2FAC*/


