
#include <R.h>
#include <Rinternals.h>

#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]

SEXP dlik (SEXP x, SEXP lx, SEXP length) {

  int i = 0, k = 0;
  unsigned int *n;
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table. */
  n = (unsigned int *) R_alloc(INT(lx), sizeof(int));
  memset(n, '\0', sizeof(int) * INT(lx));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < INT(length); k++) {

    n[INTEGER(x)[k] - 1]++;

  }/*FOR*/

  /* compute the entropy from the joint and marginal frequencies. */
  for (i = 0; i < INT(lx); i++) {

      if (n[i] != 0)
        NUM(result) += (double)n[i] * log((double)n[i] / INT(length));

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*DLIK*/

SEXP cdlik (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

  int i = 0, j = 0, k = 0;
  unsigned int **n, *nj;
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
  nj = (unsigned int *) R_alloc(INT(ly), sizeof(int));
  memset(nj, '\0', sizeof(int) * INT(ly));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < INT(length); k++) {

    n[INTEGER(x)[k] - 1][INTEGER(y)[k] - 1]++;
    nj[INTEGER(y)[k] - 1]++;

  }/*FOR*/

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < INT(lx); i++)
    for (j = 0; j < INT(ly); j++) {

      if (n[i][j] != 0)
        NUM(result) += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);

    }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CDLIK*/

