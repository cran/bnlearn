#include "common.h"

static SEXP cov2(SEXP data, SEXP length);

/* Linear Correlation, to be used in C code. */
double c_fast_cor(double *xx, double *yy, int *num) {

int i = 0;
double xm = 0, ym = 0, xsd = 0, ysd = 0, sum = 0;
double tol = MACHINE_TOL;

  /* compute the mean values.  */
  for (i = 0 ; i < *num; i++) {

    xm += xx[i];
    ym += yy[i];

  }/*FOR*/

  xm /= (*num);
  ym /= (*num);

  /* compute the actual covariance. */
  for (i = 0; i < *num; i++) {

    sum += (xx[i] - xm) * (yy[i] - ym);
    xsd += (xx[i] - xm) * (xx[i] - xm);
    ysd += (yy[i] - ym) * (yy[i] - ym);

  }/*FOR*/

  /* safety check against "divide by zero" errors. */
  if (fabs(xsd) < tol || fabs(ysd) < tol)
    sum = 0;
  else
    sum /= sqrt(xsd) * sqrt(ysd);

  /* double check that the coefficient is in the [-1, 1] range. */
  if (sum > 1) {

    warning("fixed correlation coefficient greater than 1, probably due to floating point errors.");

    sum = 1;

  }/*THEN*/
  else if (sum < -1) {

    warning("fixed correlation coefficient lesser than -1, probably due to floating point errors.");

    sum = -1;

  }/*ELSE*/

  return sum;

}/*C_FAST_COR*/

/* Linear Correlation. */
SEXP fast_cor(SEXP x, SEXP y, SEXP length) {

int *n = INTEGER(length);
double *xx = REAL(x), *yy = REAL(y);
SEXP res;

  PROTECT(res = allocVector(REALSXP, 1));
  NUM(res) = c_fast_cor(xx, yy, n);
  UNPROTECT(1);

  return res;

}/*FAST_COR*/

/* Partial Linear Correlation. */
SEXP fast_pcor(SEXP data, SEXP length, SEXP shrinkage) {

int i = 0, ncols = LENGTH(data);
int *shrink = LOGICAL(shrinkage);
double *u = NULL, *d = NULL, *vt = NULL, *res = NULL;
double k11 = 0, k12 = 0, k22 = 0;
double tol = MACHINE_TOL;
SEXP result, cov, svd;

  /* compute the covariance matrix. */
  if (*shrink > 0)
    PROTECT(cov = cov_lambda(data, length)); 
  else
    PROTECT(cov = cov2(data, length));
  /* compute the singular value decomposition of the covariance matrix. */
  PROTECT(svd = r_svd(cov));

  /* extract the three matrices form the list. */
  u = REAL(getListElement(svd, "u"));
  d = REAL(getListElement(svd, "d"));
  vt = REAL(getListElement(svd, "vt"));

  /* compute the three elements of the pseudoinverse needed
   * for the partial correlation coefficient. */
  for (i = 0; i < ncols; i++) {

    if (d[i] > tol) {

      k11 += u[CMC(0, i, ncols)] * vt[CMC(i, 0, ncols)] / d[i];
      k12 += u[CMC(0, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];
      k22 += u[CMC(1, i, ncols)] * vt[CMC(i, 1, ncols)] / d[i];

    }/*THEN*/

  }/*FOR*/

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);

  /* safety check against "divide by zero" errors. */
  if (fabs(k11) < tol || fabs(k22) < tol)
    *res = 0;
  else
    *res = -k12 / sqrt(k11 * k22);

  /* double check that partial correlation really is in [-1, 1], ill
   * conditioned matrices and numeric errors in SVD can result in
   * invalid partial correlation coefficients. */
  if (*res > 1) {

    warning("fixed partial correlation coefficient greater than 1, probably due to floating point errors.");

    *res = 1;

  }/*THEN*/
  else if (*res < -1) {

    warning("fixed partial correlation coefficient lesser than -1, probably due to floating point errors.");

    *res = -1;

  }/*ELSE*/

  UNPROTECT(3);
  return result;

}/*FAST_PCOR*/

static SEXP cov2(SEXP data, SEXP length) {

int i = 0, j = 0, k = 0, cur = 0;
int *n = INTEGER(length), ncols = LENGTH(data);
double *mean = NULL, *var = NULL, **column = NULL;
SEXP res;

  /* allocate the covariance matrix. */
  PROTECT(res = allocMatrix(REALSXP, ncols, ncols));
  var = REAL(res);
  memset(var, '\0', ncols * ncols * sizeof(double));

  /* allocate an array to store the mean values. */
  mean = Calloc(ncols, double);
  memset(mean, '\0', sizeof(double) * ncols);

  /* allocate and initialize an array of pointers for the variables. */
  column = (double **) Calloc(ncols, double *);
  for (i = 0; i < ncols; i++)
    column[i] = REAL(VECTOR_ELT(data, i));

  /* compute the mean values  */
  for (i = 0; i < ncols; i++) {

    for (j = 0 ; j < *n; j++)
      mean[i] += column[i][j];

    mean[i] /= (*n);

  }/*FOR*/

  /* compute the actual covariance. */
  for (i = 0; i < ncols; i++) {

    for (j = i; j < ncols; j++) {

      cur = CMC(i, j, ncols);

      for (k = 0; k < *n; k++)
        var[cur] += (column[i][k] - mean[i]) * (column[j][k] - mean[j]);

      var[cur] /= (*n) - 1;

      /* fill in the symmetric element of the matrix. */
      var[CMC(j, i, ncols)] = var[cur];

    }/*FOR*/

  }/*FOR*/

  Free(column);
  Free(mean);
  UNPROTECT(1);
  return res;

}/*COV2*/

