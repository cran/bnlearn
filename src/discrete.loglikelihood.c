
#include "common.h"

SEXP dlik(SEXP x, SEXP lx, SEXP length) {

int i = 0, k = 0;
int *n = NULL, *xx = INTEGER(x), *llx = INTEGER(lx), *num = INTEGER(length);
double *res = NULL;
SEXP result;

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table. */
  n = alloc1dcont(*llx);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1]++;

  }/*FOR*/

  /* compute the entropy from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++) {

    if (n[i] != 0)
      *res += (double)n[i] * log((double)n[i] / (*num));

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*DLIK*/

SEXP cdlik(SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

int i = 0, j = 0, k = 0;
int **n = NULL, *nj = NULL;
int *llx = INTEGER(lx), *lly = INTEGER(ly), *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y);
double *res = NULL;
SEXP result;

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(*llx, *lly);
  nj = alloc1dcont(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

      nj[j] += n[i][j];

    }/*FOR*/

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

      if (n[i][j] != 0)
        *res += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);

    }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CDLIK*/

