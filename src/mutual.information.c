
#include <R.h>
#include <Rinternals.h>

#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]

SEXP mi (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

  int i = 0, j = 0, k = 0;
  unsigned int n[INT(lx)][INT(ly)], ni[INT(lx)], nj[INT(ly)];
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contingency table. */
  memset(n, '\0', sizeof(int) * INT(lx) * INT(ly));

  /* initialize the marginal frequencies. */
  memset(ni, '\0', sizeof(int) * INT(lx));
  memset(nj, '\0', sizeof(int) * INT(ly));

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

SEXP cmi (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

  int i = 0, j = 0, k = 0;
  unsigned int n[INT(lx)][INT(ly)][INT(lz)], ni[INT(lx)][INT(lz)],
               nj[INT(ly)][INT(lz)], nk[INT(lz)];
  SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* initialize result to zero. */
  NUM(result) = 0;

  /* initialize the contignecy table. */
  memset(n, '\0', sizeof(int) * INT(lx) * INT(ly) * INT(lz));

  /* initialize the marginal frequencies. */
  memset(ni, '\0', sizeof(int) * INT(lx) * INT(lz));
  memset(nj, '\0', sizeof(int) * INT(ly) * INT(lz));
  memset(nk, '\0', sizeof(int) * INT(lz));

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

