#include "include/rcore.h"
#include "include/matrix.h"

/* check if a square matrix is symmetric in a zero-copy way. */
SEXP is_symmetric(SEXP matrix) {

int i = 0, j = 0, n = nrows(matrix);
double *m = NULL;

  /* dereference the matrix. */
  m = REAL(matrix);

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      if (m[CMC(i, j, n)] != m[CMC(j, i, n)])
        return ScalarLogical(FALSE);

  return ScalarLogical(TRUE);

}/*IS_SYMMETRIC*/

/* check whether a symmetric square matrix is Cauchy-Schwarz compliant. */
SEXP is_cauchy_schwarz(SEXP matrix) {

int i = 0, j = 0, n = nrows(matrix);
double *m = NULL;

  /* dereference the matrix. */
  m = REAL(matrix);

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      if (m[CMC(i, j, n)] * m[CMC(i, j, n)] > m[CMC(i, i, n)] * m[CMC(j, j, n)])
        return ScalarLogical(FALSE);

  return ScalarLogical(TRUE);

}/*IS_CAUCHY_SCHWARZ*/

