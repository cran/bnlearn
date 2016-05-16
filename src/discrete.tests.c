#include "include/rcore.h"
#include "include/tests.h"

#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP gsquare, SEXP adjusted) {

int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double *res = NULL;
SEXP result;

  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  if (isTRUE(adjusted))
    res[0] = c_chisqtest(xx, llx, yy, lly, num, res + 1, MI_ADF);
  else
    res[0] = c_chisqtest(xx, llx, yy, lly, num, res + 1, MI);

  /* rescale to match the G^2 test. */
  if (isTRUE(gsquare))
    res[0] *= 2 * num;

  UNPROTECT(1);

  return result;

}/*MI*/

/* unconditional parametric asymptotic tests for categorical data. */
double c_chisqtest(int *xx, int llx, int *yy, int lly, int num, double *df,
    test_e test) {

int  **n = NULL, *ni = NULL, *nj = NULL, adj = IS_ADF(test);
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
  fill_2d_table(xx, yy, &n, &ni, &nj, llx, lly, num);
  /* compute the mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = mi_kernel(n, ni, nj, llx, lly, num) / num;
  else if ((test == X2) || (test == X2_ADF))
    res = x2_kernel(n, ni, nj, llx, lly, num);

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? df_adjust(ni, llx, nj, lly) : (llx - 1) * (lly - 1);

  Free2D(n, llx);
  Free1D(ni);
  Free1D(nj);

  return res;

}/*C_CHISQTEST*/

/* conditional mutual information, to be used in C code. */
double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test) {

int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL, adj = IS_ADF(test);
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
  fill_3d_table(xx, yy, zz, &n, &ni, &nj, &nk, llx, lly, llz, num);
  /* compute the conditional mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = cmi_kernel(n, ni, nj, nk, llx, lly, llz) / num;
  else if ((test == X2) || (test == X2_ADF))
    res = cx2_kernel(n, ni, nj, nk, llx, lly, llz);

  /* compute the degrees of freedom. */
  if (df)
    *df = adj ? cdf_adjust(ni, llx, nj, lly, llz) : (llx - 1) * (lly - 1) * llz;

  Free3D(n, llz, llx);
  Free2D(ni, llz);
  Free2D(nj, llz);
  Free1D(nk);

  return res;

}/*C_CCHISQTEST*/

/* compute the mutual information from the joint and marginal frequencies. */
double mi_kernel(int **n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      res += MI_PART(n[i][j], nrowt[i], ncolt[j], length);

  return res;

}/*MI_KERNEL*/

/* compute Pearson's X^2 coefficient from the joint and marginal frequencies. */
double x2_kernel(int **n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length) {

int i = 0, j = 0;
double res = 0, expected = 0;

  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++) {

      expected = nrowt[i] * (double)ncolt[j] / length;

      if (expected != 0)
        res += (n[i][j] - expected) * (n[i][j] - expected) / expected;

    }/*FOR*/

  return res;

}/*X2_KERNEL*/

/* compute the conditional mutual information from the joint and marginal
 * frequencies. */
double cmi_kernel(int ***n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < nl; k++)
    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++)
        res += MI_PART(n[k][i][j], nrowt[k][i], ncolt[k][j], ncond[k]);

  return res;

}/*CMI_KERNEL*/

/* compute the Pearson's conditional X^2 coefficient from the joint and
 * marginal frequencies. */
double cx2_kernel(int ***n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl) {

int i = 0, j = 0, k = 0;
double expected = 0, res = 0;

  for (k = 0; k < nl; k++) {

    if (ncond[k] == 0) continue;

    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++) {

       expected = nrowt[k][i] * (double)ncolt[k][j] / ncond[k];

       if (expected != 0)
          res += (n[k][i][j] - expected) * (n[k][i][j] - expected) / expected;

      }/*FOR*/

  }/*FOR*/

  return res;

}/*CX2_KERNEL*/

