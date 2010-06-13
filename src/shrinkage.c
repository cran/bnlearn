
#include "common.h"

static void _mi_lambda(double *n, double *lambda, double *target, int *num,
    int *llx, int *lly, int *llz);

/* shrinked mutual information, to be used for the asymptotic test. */
SEXP shmi (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

int i = 0, j = 0, k = 0;
double **n = NULL, *ni = NULL, *nj = NULL;
int *llx = INTEGER(lx), *lly = INTEGER(ly), *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y);
double lambda = 0, target = 1/(double)((*llx) * (*lly));
double *res = NULL;
SEXP result;

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dreal(*llx, *lly);
  ni = alloc1dreal(*llx);
  nj = alloc1dreal(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* estimate the optimal lambda for the data. */
  _mi_lambda((double *)n, &lambda, &target, num, llx, lly, NULL);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
        n[i][j] = lambda * target + (1 - lambda) * n[i][j] / (*num);

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

      if (n[i][j] != 0)
        *res += n[i][j] * log(n[i][j] / (ni[i] * nj[j]));

    }/*FOR*/

  UNPROTECT(1);

  return result;

}/*SHMI*/

/* shrinked conditional mutual information, to be used for the asymptotic
 * test. */
SEXP shcmi (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

int i = 0, j = 0, k = 0;
double ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
int *llx = INTEGER(lx), *lly = INTEGER(ly), *llz = INTEGER(lz);
int *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
double lambda = 0, target = 1/(double)((*llx) * (*lly) * (*llz));
double *res = NULL;
SEXP result;

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dreal(*llx, *lly, *llz);
  ni = alloc2dreal(*llx, *llz);
  nj = alloc2dreal(*lly, *llz);
  nk = alloc1dreal(*llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;

  }/*FOR*/

  /* estimate the optimal lambda for the data. */
  _mi_lambda((double *)n, &lambda, &target, num, llx, lly, llz);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++)
        n[i][j][k] = lambda * target + (1 - lambda) * n[i][j][k] / (*num);

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  for (k = 0; k < *llz; k++) {

    /* check each level of the conditioning variable to avoid (again)
     * "divide by zero" errors. */
    if (nk[k] == 0)
      continue;

    for (j = 0; j < *lly; j++) {

      for (i = 0; i < *llx; i++) {

        if (n[i][j][k] > 0)
          *res += n[i][j][k] * log( (n[i][j][k] * nk[k]) / (ni[i][k] * nj[j][k]) );

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*SHCMI*/

/* compute the shrinkage intensity lambda for the mutual information. */
static void _mi_lambda(double *n, double *lambda, double *target, int *num,
    int *llx, int *lly, int *llz) {

double lden = 0, lnum = 0, temp = 0;

  /* compute the numerator and the denominator of the shrinkage intensity;
   * if the third dimension is a NULL pointer it's a 2-dimensional table. */
  if (!llz) {

    for (int i = 0; i < *llx; i++)
      for (int j = 0; j < *lly; j++) {

        temp = ((double **)n)[i][j] / (double)(*num);
        lnum += temp * temp;
        temp = *target - ((double **)n)[i][j] / (double)(*num);
        lden += temp * temp;

      }/*FOR*/

  }/*THEN*/
  else {

    for (int i = 0; i < *llx; i++)
      for (int j = 0; j < *lly; j++)
        for (int k = 0; k < *llz; k++) {

          temp = ((double ***)n)[i][j][k] / (double)(*num);
          lnum += temp * temp;
          temp = *target - ((double ***)n)[i][j][k] / (double)(*num);
          lden += temp * temp;

      }/*FOR*/

  }/*ELSE*/

   /* compute the shrinkage intensity (avoiding "divide by zero" errors). */
  if (lden == 0)
    *lambda = 1;
  else
    *lambda = (1 - lnum) / ((double)(*num - 1) * lden);

  /* bound the shrinkage intensity in the [0,1] interval. */
  if (*lambda > 1)
    *lambda = 1;
  if (*lambda < 0)
    *lambda = 0;

}/*_MI_LAMBDA*/


