#include "include/rcore.h"
#include "include/allocations.h"
#include "include/tests.h"

/* unconditional mutual information, to be used in C code. */
double c_mi(int *xx, int llx, int *yy, int lly, int num, double *df,
    int adj) {

int i = 0, j = 0, k = 0;
int  **n = NULL, *ni = NULL, *nj = NULL;
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

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      res += MI_PART(n[i][j], ni[i], nj[j], num);

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? df_adjust(ni, llx, nj, lly) : (llx - 1) * (lly - 1);

  return res / num;

}/*C_MI*/

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP gsquare, SEXP adjusted) {

int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double *res = NULL;
SEXP result;

  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  res[0] = c_mi(xx, llx, yy, lly, num, res + 1, isTRUE(adjusted));

  /* rescale to match the G^2 test. */
  if (isTRUE(gsquare))
    res[0] *= 2 * num;

  UNPROTECT(1);

  return result;

}/*MI*/

/* conditional mutual information, to be used in C code. */
double c_cmi(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num,
    double *df, int adj) {

int i = 0, j = 0, k = 0;
int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
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

  /* compute the conditional mutual information from the joint and
     marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++)
        res += MI_PART(n[i][j][k], ni[i][k], nj[j][k], nk[k]);

  res /= num;

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? cdf_adjust(ni, llx, nj, lly, llz) : (llx - 1) * (lly - 1) * llz;

  return res;

}/*C_CMI*/

