#include "include/rcore.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/blas.h"
#include "include/covariance.h"

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

}/*C_FAST_COR*/

/* Partial Linear Correlation. */
double c_fast_pcor(covariance cov, int v1, int v2, int *err, int decomp) {

int i = 0, errcode = 0;
double res = 0, k11 = 0, k12 = 0, k22 = 0;
double tol = MACHINE_TOL, sv_tol = 0;

  /* compute the singular value decomposition of the covariance matrix. */
  if (decomp)
    c_svd(cov.mat, cov.u, cov.d, cov.vt, &cov.dim, &cov.dim, &cov.dim,
      FALSE, &errcode);

  /* if the SVD decomposition fails, assume the partial correlation is zero. */
  if (errcode != 0) {

    if (!err)
      *err = errcode;
    else
      warning("failed to compute the pseudoinverse of the covariance matrix, assuming independence.");

    res = 0;

    return res;

  }/*THEN*/

  /* set the threshold for the singular values as in corpcor. */
  sv_tol = cov.dim * cov.d[0] * tol * tol;

  /* compute the three elements of the pseudoinverse needed
   * for the partial correlation coefficient. */
  for (i = 0; i < cov.dim; i++) {

    if (cov.d[i] > sv_tol) {

      k11 += cov.u[CMC(v1, i, cov.dim)] * cov.vt[CMC(i, v1, cov.dim)] / cov.d[i];
      k12 += cov.u[CMC(v1, i, cov.dim)] * cov.vt[CMC(i, v2, cov.dim)] / cov.d[i];
      k22 += cov.u[CMC(v2, i, cov.dim)] * cov.vt[CMC(i, v2, cov.dim)] / cov.d[i];

    }/*THEN*/

  }/*FOR*/

  /* safety check against "divide by zero" errors and negative variances. */
  res = SAFE_COR(-k12, k11, k22)
  /* double check that partial correlation is in the [-1, 1] range. */
  COR_BOUNDS(res);

  return res;

}/*C_FAST_PCOR*/

/* linear correlation from incomplete data, which also computes the number of
 * complete observations. */
double c_cor_with_missing(double *x, double *y, int nobs, double *xm,
    double *ym, double *xsd, double *ysd, int *ncomplete) {

int i = 0, nc = 0;
double tol = MACHINE_TOL;
long double xxm = 0, yym = 0, xxsd = 0, yysd = 0, cov = 0, cor = 0;

  /* compute the means. */
  for (i = 0; i < nobs; i++) {

    if (ISNAN(x[i]) || ISNAN(y[i]))
      continue;

    xxm += x[i];
    yym += y[i];
    nc++;

  }/*FOR*/

  /* if there are no complete observations, assume the correlation is zero. */
  if (nc == 0)
    goto end;

  xxm /= nc;
  yym /= nc;

  /* compute the variances and the covariance. */
  for (i = 0; i < nobs; i++) {

    if (ISNAN(x[i]) || ISNAN(y[i]))
      continue;

    xxsd += (x[i] - xxm) * (x[i] - xxm);
    yysd += (y[i] - yym) * (y[i] - yym);
    cov += (x[i] - xxm) * (y[i] - yym);

  }/*FOR*/

  /* safety check against "divide by zero" errors. */
  cor = SAFE_COR(cov, xxsd, yysd);
  /* double check that the coefficient is in the [-1, 1] range. */
  COR_BOUNDS(cor);

end:

  /* save the means, the variances and the number of complete observations. */
  *ncomplete = nc;
  if (xm)
    *xm = xxm;
  if (ym)
    *ym = yym;
  if (xsd)
    *xsd = xxsd;
  if (ysd)
    *ysd = yysd;

  return (double)cor;

}/*C_COR_WITH_MISSING*/

