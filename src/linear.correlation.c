#include "include/rcore.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/blas.h"

#define COR_BOUNDS(x) \
  if (x > 1) { \
    warning("fixed correlation coefficient greater than 1, probably due to floating point errors."); \
    x = 1; \
  }/*THEN*/ \
  else if (x < -1) { \
    warning("fixed correlation coefficient lesser than -1, probably due to floating point errors."); \
    x = -1; \
  }/*THEN*/

#define SAFE_COR(cov, xvar, yvar) \
  ((xvar < tol) || (yvar < tol)) ? 0 :  cov / sqrt(xvar * yvar);

/* linear correlation from known mean and variance, to be used in C code. */
double c_fast_cor(double *xx, double *yy, int num, double xm, double ym,
    long double xsd, long double ysd) {

int i = 0;
long double sum = 0;
double tol = MACHINE_TOL;

  /* compute the actual covariance. */
  for (i = 0; i < num; i++)
    sum += (xx[i] - xm) * (yy[i] - ym);

  /* safety check against "divide by zero" errors. */
  sum = SAFE_COR(sum, xsd, ysd);
  /* double check that the coefficient is in the [-1, 1] range. */
  COR_BOUNDS(sum);

  return (double)sum;

}/*C_FAST_COR2*/

/* Partial Linear Correlation. */
double c_fast_pcor(double *cov, double *u, double *d, double *vt, int ncols,
    int strict) {

int i = 0, errcode = 0;
double res = 0, k11 = 0, k12 = 0, k22 = 0;
double tol = MACHINE_TOL, sv_tol = 0;

  /* compute the singular value decomposition of the covariance matrix. */
  c_svd(cov, u, d, vt, &ncols, &ncols, &ncols, FALSE, &errcode);

  if (errcode != 0) {

    if (strict) {

      error("failed to compute the pseudoinverse of the covariance matrix.");

    }/*THEN*/
    else {

      /* if computing SVD decomposition fails, assume a null correlation. */
      res = 0;
      /* warn the user that something went wrong.*/
      warning("failed to compute the pseudoinverse of the covariance matrix, assuming independence.");

      return res;

    }/*ELSE*/

  }/*THEN*/

  /* set the threshold for the singular values as in corpcor. */
  sv_tol = ncols * d[0] * tol * tol;

  /* compute the three elements of the pseudoinverse needed
   * for the partial correlation coefficient. */
  for (i = 0; i < ncols; i++) {

    if (d[i] > sv_tol) {

      k11 += u[CMC(0, i, ncols)] * vt[CMC(i, 0, ncols)] / d[i];
      k12 += u[CMC(0, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];
      k22 += u[CMC(1, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];

    }/*THEN*/

  }/*FOR*/

  /* safety check against "divide by zero" errors and negative variances. */
  res = SAFE_COR(-k12, k11, k22)
  /* double check that partial correlation is in the [-1, 1] range. */
  COR_BOUNDS(res);

  return res;

}/*C_FAST_PCOR*/

