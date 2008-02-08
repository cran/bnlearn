
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>

void dlik (int *x, int *lx, int *length, double *result) {

  int i = 0, k = 0;
  unsigned int n[*lx];

  /* initialize result to zero. */
  *result = 0;

  /* initialize the contingency table. */
  memset(n, '\0', sizeof(int)*(*lx));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *length; k++) {

    n[x[k] - 1]++;

  }/*FOR*/

  for (i = 0; i < *lx; i++) {

      if (n[i] != 0)
        *result += (double)n[i] * log((double)n[i] / (double)*length);

  }/*FOR*/

}/*DLIK*/

void cdlik (int *x, int *y, int *lx, int *ly, int *length, double *result) {

  int i = 0, j = 0, k = 0;
  unsigned int n[*lx][*ly], nj[*ly];

  /* initialize result to zero. */
  *result = 0;

  /* initialize the contingency table. */
  memset(n, '\0', sizeof(int)*(*lx)*(*ly));

  /* initialize the marginal frequencies. */
  memset(nj, '\0', sizeof(int)*(*ly));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *length; k++) {

    n[x[k] - 1][y[k] - 1]++;
    nj[y[k] - 1]++;

  }/*FOR*/

  for (i = 0; i < *lx; i++)
    for (j = 0; j < *ly; j++) {

      if (n[i][j] != 0)
        *result += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);

    }/*FOR*/

}/*CDLIK*/

