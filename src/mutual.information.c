
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

void mi (int *x, int *y, int *lx, int *ly, int *length, double *result) {

  int i = 0, j = 0, k = 0;
  unsigned int n[*lx][*ly], ni[*lx], nj[*ly];

  /* initialize result to zero. */
  *result = 0;

  /* initialize the contignecy table. */
  memset(n, '\0', sizeof(int)*(*lx)*(*ly));

  /* initialize the marginal frequencies. */
  memset(ni, '\0', sizeof(int)*(*lx));
  memset(nj, '\0', sizeof(int)*(*ly));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *length; k++) {

    n[x[k] - 1][y[k] - 1]++;
    ni[x[k] - 1]++;
    nj[y[k] - 1]++;

  }/*FOR*/

  for (i = 0; i < *lx; i++)
    for (j = 0; j < *ly; j++) {

      if (n[i][j] != 0)
        *result += ((double)n[i][j]) *
                log((double)n[i][j]*(*length)/(double)(ni[i]*nj[j]));

    }/*FOR*/

    *result = *result/(*length);

}/*MI*/

void cmi (int *x, int *y, int *z, int *lx, int *ly, int *lz, int *length, double *result) {

  int i = 0, j = 0, k = 0;
  unsigned int n[*lx][*ly][*lz], ni[*lx][*lz], nj[*ly][*lz], nk[*lz];

  /* initialize result to zero. */
  *result = 0;

  /* initialize the contignecy table. */
  memset(n, '\0', sizeof(int)*(*lx)*(*ly)*(*lz));

  /* initialize the marginal frequencies. */
  memset(ni, '\0', sizeof(int)*(*lx)*(*lz));
  memset(nj, '\0', sizeof(int)*(*ly)*(*lz));
  memset(nk, '\0', sizeof(int)*(*lz));

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < *length; k++) {

    n[x[k] - 1][y[k] - 1][z[k] - 1]++;
    ni[x[k] - 1][z[k] - 1]++;
    nj[y[k] - 1][z[k] - 1]++;
    nk[z[k] - 1]++;

  }/*FOR*/

  for (i = 0; i < *lx; i++)
    for (j = 0; j < *ly; j++)
      for (k = 0; k < *lz; k++) {

       if (n[i][j][k] != 0) {

          *result += (double)n[i][j][k] *
            log( (double)(n[i][j][k]*nk[k]) / (double)(ni[i][k]*nj[j][k]) );

        }/*THEN*/

      }/*FOR*/

   *result = *result/(*length);

}/*CMI*/

