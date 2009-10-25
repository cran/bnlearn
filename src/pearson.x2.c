
#include "common.h"

SEXP x2 (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

int i = 0, j = 0, k = 0, **n = NULL, *ni = NULL, *nj = NULL;
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
  ni = alloc1dcont(*llx);
  nj = alloc1dcont(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

      ni[i] += n[i][j];
      nj[j] += n[i][j];

    }/*FOR*/

  /* compute the X^2 from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

      if (n[i][j] != 0)
        *res += (n[i][j] - ni[i] * (double)nj[j] / (*num)) *
                (n[i][j] - ni[i] * (double)nj[j] / (*num)) /
                (ni[i] * (double)nj[j] / (*num));

    }/*FOR*/

  UNPROTECT(1);

  return result;

}/*X2*/

SEXP cx2 (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

int i = 0, j = 0, k = 0, ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
int *llx = INTEGER(lx), *lly = INTEGER(ly), *llz = INTEGER(lz), *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
double *res = NULL;
SEXP result;

  /* allocate  and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dcont(*llx, *lly, *llz);
  ni = alloc2dcont(*llx, *llz);
  nj = alloc2dcont(*lly, *llz);
  nk = alloc1dcont(*llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  /* compute the conditional X^2 from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) {

        if (n[i][j][k] != 0)
          *res += (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) *
                  (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) /
                  (ni[i][k] * (double)nj[j][k] / nk[k]);

      }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CX2*/

