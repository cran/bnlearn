
#include <R.h>
#include <Rinternals.h>

#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]

SEXP x2 (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

  int i = 0, j = 0, k = 0;
  unsigned int **n, *ni, * nj;
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table. */
  n = (unsigned int **) R_alloc(INT(lx), sizeof(int *));
  for (i = 0; i < INT(lx); i++) {

    n[i] = (unsigned int *) R_alloc(INT(ly), sizeof(int));
    memset(n[i], '\0', sizeof(int) * INT(ly));

  }/*FOR*/

  /* initialize the marginal frequencies. */
  ni = (unsigned int *) R_alloc(INT(lx), sizeof(int));
  memset(ni, '\0', sizeof(int) * INT(lx));
  nj = (unsigned int *) R_alloc(INT(ly), sizeof(int));
  memset(nj, '\0', sizeof(int) * INT(ly));

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
  unsigned int ***n, **ni, **nj, *nk;
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table. */
  n = (unsigned int ***) R_alloc(INT(lx), sizeof(int *));
  for (i = 0; i < INT(lx); i++) {

    n[i] = (unsigned int **) R_alloc(INT(ly), sizeof(int *));

    for (j = 0; j < INT(ly); j++) {

      n[i][j] = (unsigned int*) R_alloc(INT(lz), sizeof(int));
      memset(n[i][j], '\0', sizeof(int) * INT(lz));

    }/*FOR*/

  }/*FOR*/

  /* initialize the marginal frequencies. */
  ni = (unsigned int **) R_alloc(INT(lx), sizeof(int *));
  for (i = 0; i < INT(lx); i++) {

    ni[i] = (unsigned int *) R_alloc(INT(lz), sizeof(int));
    memset(ni[i], '\0', sizeof(int) * INT(lz));

  }/*FOR*/

  nj = (unsigned int **) R_alloc(INT(ly), sizeof(int *));
  for (i = 0; i < INT(ly); i++) {

    nj[i] = (unsigned int *) R_alloc(INT(lz), sizeof(int));
    memset(nj[i], '\0', sizeof(int) * INT(lz));

  }/*FOR*/

  nk = (unsigned int *) R_alloc(INT(lz), sizeof(int));
  memset(nk, '\0', sizeof(int) * INT(lz));

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

