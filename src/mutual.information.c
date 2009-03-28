
#include "common.h"

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

  int i = 0, j = 0, k = 0;
  int **n, *ni, *nj;
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

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < INT(lx); i++)
    for (j = 0; j < INT(ly); j++) {

      if (n[i][j] != 0)
        NUM(result) += ((double)n[i][j]) *
                log((double)n[i][j]*(INT(length))/(double)(ni[i]*nj[j]));

    }/*FOR*/

  NUM(result) = NUM(result)/INT(length);

  UNPROTECT(1);

  return result;

}/*MI*/

/* conditional mutual information, to be used for the asymptotic test. */
SEXP cmi (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

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

  /* compute the conditional mutual information from the joint and
     marginal frequencies. */
  for (i = 0; i < INT(lx); i++)
    for (j = 0; j < INT(ly); j++)
      for (k = 0; k < INT(lz); k++) {

       if (n[i][j][k] != 0) {

          NUM(result) += (double)n[i][j][k] *
            log( (double)(n[i][j][k]*nk[k]) / (double)(ni[i][k]*nj[j][k]) );

        }/*THEN*/

      }/*FOR*/

  NUM(result) = NUM(result)/INT(length);

  UNPROTECT(1);

  return result;

}/*CMI*/

