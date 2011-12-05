#include "common.h"

/* check if a square matrix is symmetric in a zero-copy way. */
SEXP is_symmetric(SEXP matrix) {

int i = 0, j = 0, n = nrows(matrix);
double *m = NULL;
SEXP result;

  /* dereference the matrix. */
  m = REAL(matrix);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(LGLSXP, 1));
  LOGICAL(result)[0] = TRUE;

  for (i = 0; i < n; i++) {

    for (j = i + 1; j < n; j++) {

      /* two cells do no match; it's useless to go on, set the return
       * value to FALSE and jump at the end of the function. */
      if (m[CMC(i, j, n)] != m[CMC(j, i, n)]) {

        LOGICAL(result)[0] = FALSE;

        goto end;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

end:

  UNPROTECT(1);

  return result;

}/*IS_SYMMETRIC*/

/* check whether a symmetric square matrix is Cauchy-Schwarz compliant. */
SEXP is_cauchy_schwarz(SEXP matrix) {

int i = 0, j = 0, n = nrows(matrix);
double *m = NULL;
SEXP result;

  /* dereference the matrix. */
  m = REAL(matrix);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(LGLSXP, 1));
  LOGICAL(result)[0] = TRUE;

  for (i = 0; i < n; i++) {

    for (j = i + 1; j < n; j++) {

      if (m[CMC(i, j, n)] * m[CMC(i, j, n)] > m[CMC(i, i, n)] * m[CMC(j, j, n)]) {

        LOGICAL(result)[0] = FALSE;

        goto end;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

end:

  UNPROTECT(1);

  return result;

}/*IS_CAUCHY_SCHWARZ*/

