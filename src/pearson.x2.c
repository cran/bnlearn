
#include "common.h"

SEXP x2 (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

  int i = 0, j = 0, k = 0;
  int **n, *ni, * nj;
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(INT(lx), INT(ly));
  ni = alloc1dcont(INT(lx));
  nj = alloc1dcont(INT(ly));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < INT(length); k++) {

    n[INTEGER(x)[k] - 1][INTEGER(y)[k] - 1]++;
    ni[INTEGER(x)[k] - 1]++;
    nj[INTEGER(y)[k] - 1]++;

  }/*FOR*/

  /* compute the X^2 from the joint and marginal frequencies. */
  for (i = 0; i < INT(lx); i++)
    for (j = 0; j < INT(ly); j++) {

      if (n[i][j] != 0)
        NUM(result) += (n[i][j] - ni[i] * (double)nj[j] / INT(length)) *
                       (n[i][j] - ni[i] * (double)nj[j] / INT(length)) /
                       (ni[i] * (double)nj[j] / INT(length));

    }/*FOR*/

  UNPROTECT(1);

  return result;

}/*X2*/

SEXP cx2 (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

  int i = 0, j = 0, k = 0;
  int ***n, **ni, **nj, *nk;
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dcont(INT(lx), INT(ly), INT(lz));
  ni = alloc2dcont(INT(lx), INT(lz));
  nj = alloc2dcont(INT(ly), INT(lz));
  nk = alloc1dcont(INT(lz));

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < INT(length); k++) {

    n[INTEGER(x)[k] - 1][INTEGER(y)[k] - 1][INTEGER(z)[k] - 1]++;
    ni[INTEGER(x)[k] - 1][INTEGER(z)[k] - 1]++;
    nj[INTEGER(y)[k] - 1][INTEGER(z)[k] - 1]++;
    nk[INTEGER(z)[k] - 1]++;

  }/*FOR*/

  /* compute the conditional X^2 from the joint and marginal frequencies. */
  for (i = 0; i < INT(lx); i++)
    for (j = 0; j < INT(ly); j++)
      for (k = 0; k < INT(lz); k++) {

        if (n[i][j][k] != 0)
          NUM(result) += (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) *
                         (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) /
                         (ni[i][k] * (double)nj[j][k] / nk[k]);

      }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CX2*/

