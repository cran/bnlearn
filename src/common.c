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

  int i, j;

  for (i = 0; i < n; i++)
    x[i] = i;
  for (i = 0; i < k; i++) {

    j = n * unif_rand();
    y[i] = x[j] + 1;
    x[j] = x[--n];

  }/*FOR*/

}/*SAMPLENOREPLACE*/

