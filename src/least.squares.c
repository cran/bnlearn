#include "include/rcore.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/globals.h"
#include "include/matrix.h"

/* fast implementation of least squares with just the intercept. */
static void c_ols0(double *y, int nrow, double *fitted, double *resid,
    double *beta, double *sd) {

int i = 0;
double mean = c_mean(y, nrow);

   /* the only coefficient is the intercept, which is the mean of the response. */
  if (beta)
    *beta = mean;
  if (sd)
    c_sd(y, nrow, 1, mean, FALSE, sd);
  if (fitted)
    for (i = 0; i < nrow; i++)
      fitted[i] = mean;
  if (resid)
    for (i  = 0; i < nrow; i++)
      resid[i] = y[i] - mean;

}/*C_OLS0*/

/* fast implementation of least squares with one regression coefficient. */
static void c_ols1(double **x, double *y, int nrow, double *fitted, double *resid,
    double *beta, double *sd) {

int i = 0, singular = FALSE;
double *m[2] = {y, *x}, mean[2] = {0, 0}, cov[4] = {0, 0, 0, 0}, a = 0, b = 0;
long double sse = 0;

  /* the regression coefficients are computed using the closed form estimators
   * for simple regression. */
  c_meanvec(m, mean, nrow, 2, 0);
  c_covmat(m, mean, nrow, 2, cov, 0);

  /* the only way for a simple regression to be singular is if the regressor
   * is constant, so it is collinear with the intercept. */
  singular = (fabs(cov[3]) < MACHINE_TOL);

  if (singular) {

    b = NA_REAL;
    a = mean[0];

  }/*THEN*/
  else {

    b = cov[1] / cov[3];
    a = mean[0] - mean[1] * b;

  }/*ELSE*/

  if (beta) {

    beta[0] = a;
    beta[1] = b;

  }/*THEN*/

  if (fitted) {

    if (singular)
      for (i = 0; i < nrow; i++)
        fitted[i] = a;
    else
      for (i = 0; i < nrow; i++)
        fitted[i] = a + b * x[0][i];

  }/*THEN*/

  if (resid) {

    if (fitted)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - fitted[i];
    else {

      if (singular)
        for (i  = 0; i < nrow; i++)
          resid[i] = y[i] - a + b * x[0][i];
      else
        for (i  = 0; i < nrow; i++)
          resid[i] = y[i] - a;

    }/*ELSE*/

  }/*THEN*/

  if (sd) {

    SD_GUARD(nrow, 2, *sd,
      if (resid)
        for (i = 0; i < nrow; i++)
          sse += resid[i] * resid[i];
      else if (fitted)
        for (i = 0; i < nrow; i++)
          sse += (y[i] - fitted[i]) * (y[i] - fitted[i]);
      else {

        if (singular)
          sse = c_sse(y, a, nrow);
        else
          for (i = 0; i < nrow; i++)
            sse += (y[i] - a - b * x[0][i]) * (y[i] - a - b * x[0][i]);

      }/*ELSE*/

      *sd = sqrt(sse / (nrow - 2));
    )

  }/*THEN*/

}/*C_OLS1*/

/* compute least squares efficiently by special-casing whenever possible. */
void c_ols(double **x, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd) {

double *qr = NULL;

  if (ncol == 0) {

    /* special case: null model. */
    c_ols0(y, nrow, fitted, resid, beta, sd);

  }/*THEN*/
  else if (ncol == 1) {

    /* special case: simple regression. */
    c_ols1(x, y, nrow, fitted, resid, beta, sd);

  }/*THEN*/
  else {

    /* prepare the data for the QR decomposition. */
    qr = Calloc1D(nrow * (ncol + 1), sizeof(double));
    c_qr_matrix(qr, x, nrow, ncol);
    /* general case: multiple regression. */
    c_qr(qr, y, nrow, ncol + 1, fitted, resid, beta, sd);

    Free1D(qr);

  }/*ELSE*/

}/*C_OLS*/

/* fast implementation of conditional least squares with just the intercept. */
static void c_cls0(double *y, int *z, int nrow, int ncond, double *fitted,
    double *resid, double *beta, double *sd) {

int i = 0, *nz = NULL;
long double *means = NULL;

  /* the matrix of regression coefficients contains the conditional mean given
   * each configuration; estimate in a single pass (set to NaN if there a no
   * observations for a particular observation). */
  means = Calloc1D(ncond, sizeof(long double));
  nz = Calloc1D(ncond, sizeof(int));

  for (i = 0; i < nrow; i++) {

    means[z[i] - 1] += y[i];
    nz[z[i] - 1]++;

  }/*FOR*/

  for (i = 0; i < ncond; i++) {

    if (nz[i] == 0)
      means[i] = R_NaN;
    else
      means[i] /= nz[i];

  }/*FOR*/

  if (beta)
    for (i  = 0; i < ncond; i++)
      beta[i] = means[i];
  if (fitted)
    for (i  = 0; i < nrow; i++)
      fitted[i] = means[z[i] - 1];
  if (resid)
    for (i  = 0; i < nrow; i++)
      resid[i] = y[i] - means[z[i] - 1];
  if (sd)
    c_cgsd(y, z, nz, nrow, ncond, 1, means, sd);

  Free1D(nz);
  Free1D(means);

}/*C_CLS0*/

/* fast implementation of conditional least squares with one regression
 * coefficient. */
static void c_cls1(double **x, double *y, int *z, int nrow, int ncond,
    double *fitted, double *resid, double *beta, double *sd) {

int i = 0, *nz = NULL;
double *cc = NULL;
long double *mean_y = NULL, *mean_x = NULL, *var_x = NULL, *cov = NULL, *ssr = NULL;

  /* first pass: compute the conditional means of y and x. */
  mean_y = Calloc1D(ncond, sizeof(long double));
  mean_x = Calloc1D(ncond, sizeof(long double));
  var_x = Calloc1D(ncond, sizeof(long double));
  cov = Calloc1D(ncond, sizeof(long double));
  nz = Calloc1D(ncond, sizeof(int));
  cc = Calloc1D(2 * ncond, sizeof(double));

  for (i = 0; i < nrow; i++) {

    mean_y[z[i] - 1] += y[i];
    mean_x[z[i] - 1] += x[0][i];
    nz[z[i] - 1]++;

  }/*FOR*/

  for (i = 0; i < ncond; i++) {

    mean_y[i] /= nz[i];
    mean_x[i] /= nz[i];

  }/*FOR*/

  /* second pass: compute the covariance between x and y and the variance of x. */
  for (i = 0; i < nrow; i++) {

    var_x[z[i] - 1] += (x[0][i] - mean_x[z[i] - 1]) * (x[0][i] - mean_x[z[i] - 1]);
    cov[z[i] - 1] += (y[i] - mean_y[z[i] - 1]) * (x[0][i] - mean_x[z[i] - 1]);

  }/*FOR*/

  /* the regression coefficients are computed using the closed form estimators
   * for simple regression (set to NaN if there a no observations for a
   * particular observation). */
  for (i = 0; i < ncond; i++) {

    if (nz[i] == 0) {

      cc[i * 2] = cc[i * 2 + 1] = R_NaN;

    }/*THEN*/
    else {

      cc[i * 2 + 1] = (fabs(var_x[i]) < MACHINE_TOL) ? 0 : cov[i] / var_x[i];
      cc[i * 2] = mean_y[i] - mean_x[i] * cc[i * 2 + 1];

    }/*ELSE*/

  }/*FOR*/

  if (beta)
    memcpy(beta, cc, 2 * ncond * sizeof(double));
  if (fitted)
    for (i  = 0; i < nrow; i++)
      fitted[i] = cc[(z[i] - 1) * 2] + x[0][i] * cc[(z[i] - 1) * 2 + 1];
  if (resid) {

    if (fitted)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - fitted[i];
    else
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - cc[(z[i] - 1) * 2] - x[0][i] * cc[(z[i] - 1) * 2 + 1];

  }/*THEN*/
  if (sd) {

    ssr = Calloc1D(ncond, sizeof(long double));

    if (resid)
      for (i  = 0; i < nrow; i++)
        ssr[z[i] - 1] += resid[i] * resid[i];
    else if (fitted)
      for (i  = 0; i < nrow; i++)
        ssr[z[i] - 1] += (y[i] - fitted[i]) * (y[i] - fitted[i]);
    else
      for (i  = 0; i < nrow; i++)
        ssr[z[i] - 1] += (y[i] - cc[(z[i] - 1) * 2] - x[0][i] * cc[(z[i] - 1) * 2 + 1]) *
               (y[i] - cc[(z[i] - 1) * 2] - x[0][i] * cc[(z[i] - 1) * 2 + 1]);

    for (i = 0; i < ncond; i++)
      SD_GUARD(nz[i], 2, sd[i],
        sd[i] = sqrt(ssr[i] / (nz[i] - 2));
      )

    Free1D(ssr);

  }/*THEN*/

  Free1D(cc);
  Free1D(mean_y);
  Free1D(mean_x);
  Free1D(var_x);
  Free1D(cov);
  Free1D(nz);

}/*C_CLS1*/

/* general implementation of conditional least squares. */
void c_clsp(double **x, double *y, int *z, int nrow, int ncol, int ncond,
    double *fitted, double *resid, double *beta, double *sd) {

int i = 0, j = 0, *nz = NULL, *counter = NULL, qr_ncol = ncol + 1, **zid = NULL;
double **qr_matrix = NULL, **yy = NULL, *ff = NULL, *rr = NULL;

  /* first pass: count how many observations for each conditioning. */
  nz = Calloc1D(ncond, sizeof(int));
  for (i = 0; i < nrow; i++)
    nz[z[i] - 1]++;

  /* allocate the responses and the matrices to pass to c_qr(). */
  yy = Calloc1D(ncond, sizeof(double *));
  for (i = 0; i < ncond; i++)
    yy[i] = Calloc1D(nz[i], sizeof(double));

  qr_matrix = Calloc1D(ncond, sizeof(double *));
  for (i = 0; i < ncond; i++)
    qr_matrix[i] = Calloc1D(nz[i] * qr_ncol, sizeof(double));

  /* set the intercept. */
  for (i = 0; i < ncond; i++)
   for (j = 0; j < nz[i]; j++)
     qr_matrix[i][j] = 1;

  /* second pass: copy the data into the matrices. */
  counter = Calloc1D(ncond, sizeof(int));

  for (i = 0; i < nrow; i++) {

    for (j = 1; j < qr_ncol; j++)
      qr_matrix[z[i] - 1][CMC(counter[z[i] - 1], j, nz[z[i] - 1])] = x[j - 1][i];

    yy[z[i] - 1][counter[z[i] - 1]] = y[i];

    counter[z[i] - 1]++;

  }/*FOR*/

  /* track the indexes of the observations in each conditioning, if needed. */
  if (fitted || resid) {

    zid = Calloc1D(ncond, sizeof(int *));
    for (i = 0; i < ncond; i++)
      zid[i] = Calloc1D(nz[i], sizeof(int));

    memset(counter, '\0', ncond * sizeof(int));
    for (i = 0; i < nrow; i++)
      zid[z[i] - 1][counter[z[i] - 1]++] = i;

  }/*THEN*/

  /* iterate over each conditioning and compute the QR decomposition. */
  memset(counter, '\0', ncond * sizeof(int));

  for (i = 0; i < ncond; i++) {

    /* special case: unobserved conditioning means everything in NaN. */
    if (nz[i] == 0) {

      if (beta)
        for (j = 0; j < qr_ncol; j++)
          beta[i * qr_ncol + j] = R_NaN;
      if (sd)
        sd[i] = R_NaN;
      continue;

    }/*THEN*/

    /* gcc cannot optimize the copying of fitted values and residuals if
     * implemented as a loop with lots of conditionals, so treat each case
     * separately. */
    if (!fitted && !resid) {

      c_qr(qr_matrix[i], yy[i], nz[i], qr_ncol, NULL, NULL,
        (!beta) ? NULL : beta + i * qr_ncol, (!sd) ? NULL : sd + i);

      continue;

    }/*THEN*/
    else if (fitted && !resid) {

      ff = Calloc1D(nz[i], sizeof(double));

      c_qr(qr_matrix[i], yy[i], nz[i], qr_ncol, ff, NULL,
        (!beta) ? NULL : beta + i * qr_ncol, (!sd) ? NULL : sd + i);

      for (j = 0; j < nz[i]; j++)
        fitted[zid[i][j]] = ff[j];

      Free1D(ff);

    }/*THEN*/
    else if (!fitted && resid) {

      rr = Calloc1D(nz[i], sizeof(double));

      c_qr(qr_matrix[i], yy[i], nz[i], qr_ncol, NULL, rr,
        (!beta) ? NULL : beta + i * qr_ncol, (!sd) ? NULL : sd + i);

      for (j = 0; j < nz[i]; j++)
        resid[zid[i][j]] = rr[j];

      Free1D(rr);

    }/*THEN*/
    else {

      ff = Calloc1D(nz[i], sizeof(double));
      rr = Calloc1D(nz[i], sizeof(double));

      c_qr(qr_matrix[i], yy[i], nz[i], qr_ncol, ff, rr,
        (!beta) ? NULL : beta + i * qr_ncol, (!sd) ? NULL : sd + i);

      for (j = 0; j < nz[i]; j++) {

        fitted[zid[i][j]] = ff[j];
        resid[zid[i][j]] = rr[j];

      }/*FOR*/

      Free1D(ff);
      Free1D(rr);

    }/*ELSE*/

  }/*FOR*/

  /* post-process the regression coefficients and replace NAs with zeros, while
   * leaving NaNs in place for unobserved conditionings. */
  if (beta)
    for (i = 0; i < ncond * qr_ncol; i++)
      if (ISNA(beta[i]))
        beta[i] = 0;

  Free1D(nz);
  Free1D(counter);
  Free2D(yy, ncond);
  Free2D(qr_matrix, ncond);
  if (fitted || resid)
    Free2D(zid, ncond);

}/*C_CLSP*/

/* compute conditional least squares efficiently by special-casing whenever
 * possible. */
void c_cls(double **x, double *y, int *z, int nrow, int ncol, int ncond,
    double *fitted, double *resid, double *beta, double *sd) {

  if (ncol == 0) {

    c_cls0(y, z, nrow, ncond, fitted, resid, beta, sd);

  }/*THEN*/
  else if (ncol == 1) {

   c_cls1(x, y, z, nrow, ncond, fitted, resid, beta, sd);

  }/*THEN*/
  else {

    c_clsp(x, y, z, nrow, ncol, ncond, fitted, resid, beta, sd);

  }/*ELSE*/

}/*C_CLS*/


