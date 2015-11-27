#include "include/rcore.h"

/* adjust degrees of freedom: zeroes are considered structural if they
 * are part of a column or a row with a zero marginal, the rest are
 * considered sampling zeros. */
double df_adjust(int *ni, int llx, int *nj, int lly) {

int i = 0, j = 0, alx = 0, aly = 0;

  for (i = 0; i < llx; i++)
    alx += (ni[i] > 0);
  for (j = 0; j < lly; j++)
    aly += (nj[j] > 0);

  /* ensure the degrees of freedom will not be negative. */
  alx = (alx >= 1) ? alx : 1;
  aly = (aly >= 1) ? aly : 1;

  return (double)((alx - 1) * (aly - 1));

}/*DF_ADJUST*/

/* adjust degrees of freedom: same as the unconditional case, but stacked
 * across the strata of the conditioning variables. */
double cdf_adjust(int **ni, int llx, int **nj, int lly, int llz) {

int i = 0, j = 0, k = 0, alx = 0, aly = 0;
double df = 0;

  for (k = 0; k < llz; k++) {

    alx = aly = 0;

    for (i = 0; i < llx; i++)
      alx += (ni[k][i] > 0);
    for (j = 0; j < lly; j++)
      aly += (nj[k][j] > 0);

    /* ensure the degrees of freedom will not be negative. */
    alx = (alx >= 1) ? alx : 1;
    aly = (aly >= 1) ? aly : 1;

    df += (alx - 1) * (aly - 1);

  }/*FOR*/

  return df;

}/*CDF_ADJUST*/

