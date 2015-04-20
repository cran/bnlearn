#include "include/rcore.h"
#include "include/allocations.h"
#include "include/tests.h"

double c_x2(int *xx, int llx, int *yy, int lly, int num, double *df,
    int adj) {

int i = 0, j = 0, k = 0, **n = NULL, *ni = NULL, *nj = NULL;
double res = 0;

  if (adj) {

    /* if there are less than 5 observations per cell on average, assume the
     * test does not have enough power and return independence. */
    if (num < 5 * llx * lly) {

      if (df) *df = 1;

      return 0;

    }/*THEN*/

  }/*THEN*/

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(llx, lly);
  ni = alloc1dcont(llx);
  nj = alloc1dcont(lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

      ni[i] += n[i][j];
      nj[j] += n[i][j];

    }/*FOR*/

  /* compute the X^2 from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

      if (n[i][j] != 0)
        res += (n[i][j] - ni[i] * (double)nj[j] / num) *
               (n[i][j] - ni[i] * (double)nj[j] / num) /
               (ni[i] * (double)nj[j] / num);

    }/*FOR*/

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? df_adjust(ni, llx, nj, lly) : (llx - 1) * (lly - 1);

  return res;

}/*C_X2*/

double c_cx2(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num,
    double *df, int adj) {

int i = 0, j = 0, k = 0, ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
double res = 0;

  if (adj) {

    /* if there are less than 5 observations per cell on average, asuume the
     * test does not have enough power and return independence. */
    if (num < 5 * llx * lly * llz) {

      if (df) *df = 1;

      return 0;

    }/*THEN*/

  }/*THEN*/

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dcont(llx, lly, llz);
  ni = alloc2dcont(llx, llz);
  nj = alloc2dcont(lly, llz);
  nk = alloc1dcont(llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  /* compute the conditional X^2 from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        if (n[i][j][k] != 0)
          res += (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) *
                 (n[i][j][k] - ni[i][k] * (double)nj[j][k] / nk[k]) /
                 (ni[i][k] * (double)nj[j][k] / nk[k]);

      }/*FOR*/

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? cdf_adjust(ni, llx, nj, lly, llz) : (llx - 1) * (lly - 1) * llz;

  return res;

}/*C_CX2*/

