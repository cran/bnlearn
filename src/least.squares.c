#include "include/rcore.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/globals.h"

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
double *m[2] = {y, *x}, mean[2] = {0, 0}, a = 0, b = 0;
long double sse = 0;
covariance cov = { .dim = 2, .mat = (double[]){ 0, 0, 0, 0 },
                   .u = NULL, .d = NULL, .vt = NULL };

  /* the regression coefficients are computed using the closed form estimators
   * for simple regression. */
  c_meanvec(m, mean, nrow, 2, 0);
  c_covmat(m, mean, nrow, 2, cov, 0);

  /* the only way for a simple regression to be singular is if the regressor
   * is constant, so it is collinear with the intercept. */
  singular = (fabs(cov.mat[3]) < MACHINE_TOL);

  if (singular) {

    b = NA_REAL;
    a = mean[0];

  }/*THEN*/
  else {

    b = cov.mat[1] / cov.mat[3];
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
          resid[i] = y[i] - a - b * x[0][i];
      else
        for (i  = 0; i < nrow; i++)
          resid[i] = y[i] - a;

    }/*ELSE*/

  }/*THEN*/

  if (sd) {

    if (nrow == 0)
      *sd = R_NaN;
    else if (nrow <= 2)
      *sd = 0;
    else {

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

    }/*ELSE*/

  }/*THEN*/

}/*C_OLS1*/

/* fast implementation of least squares with two regression coefficients. */
static void c_ols2(double **x, double *y, int nrow, double *fitted, double *resid,
    double *beta, double *sd) {

int i = 0, singular1 = FALSE, singular2 = FALSE;
double *m[3] = {y, *x, *(x + 1)}, a = 0, b1 = 0, b2 = 0, den = 0;
double mean[3] = {0, 0, 0};
long double sse = 0;
covariance cov = { .dim = 3, .mat = (double[]){ 0, 0, 0, 0, 0, 0, 0, 0, 0 },
                   .u = NULL, .d = NULL, .vt = NULL };

  /* the regression coefficients are computed using the closed form estimates
   * for the two-variable regression, which are effectively the same as the
   * corresponding partial correlation estimates. */
  c_meanvec(m, mean, nrow, 3, 0);
  c_covmat(m, mean, nrow, 3, cov, 0);

  /* there are three possible collinear configurations:
   *   1) the first variable is constant and collinear with the response;
   *   2) the second variable is constant and collinear with the response;
   *   3) the two variables are collinear with each other. */
  singular1 = (fabs(cov.mat[4]) < MACHINE_TOL);
  singular2 = (fabs(cov.mat[8]) < MACHINE_TOL) ||
              (fabs(cov.mat[5]) / sqrt(cov.mat[4] * cov.mat[8]) > 1 - MACHINE_TOL);

  if (singular1 && !singular2) {

    b1 = NA_REAL;
    b2 = cov.mat[2] / cov.mat[8];
    a =  mean[0] - b2 * mean[2];

  }/*THEN*/
  else if (!singular1 && singular2) {

    b1 = cov.mat[1] / cov.mat[4];
    b2 = NA_REAL;
    a =  mean[0] - b1 * mean[1];

  }/*THEN*/
  else if (singular1 && singular2) {

    b1 = b2 = NA_REAL;
    a = mean[0];

  }/*THEN*/
  else {

    den = (cov.mat[4] * cov.mat[8] - cov.mat[5] * cov.mat[5]);
    b1 = (cov.mat[8] * cov.mat[1] - cov.mat[5] * cov.mat[2]) / den;
    b2 = (cov.mat[4] * cov.mat[2] - cov.mat[5] * cov.mat[1]) / den;
    a = mean[0] - b1 * mean[1] - b2 * mean[2];

  }/*THEN*/

  if (beta) {

    beta[0] = a;
    beta[1] = b1;
    beta[2] = b2;

  }/*THEN*/

  if (fitted) {

    if (singular1 && !singular2)
      for (i = 0; i < nrow; i++)
        fitted[i] = a + b2 * x[1][i];
    else if (!singular1 && singular2)
      for (i = 0; i < nrow; i++)
        fitted[i] = a + b1 * x[0][i];
    else if (singular1 && singular2)
      for (i = 0; i < nrow; i++)
        fitted[i] = a;
    else
      for (i = 0; i < nrow; i++)
        fitted[i] = a + b1 * x[0][i] + b2 * x[1][i];

  }/*THEN*/

  if (resid) {

    if (fitted)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - fitted[i];
    else {

    if (singular1 && !singular2)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - a - b2 * x[1][i];
    else if (!singular1 && singular2)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - a - b1 * x[0][i];
    else if (singular1 && singular2)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - a;
    else
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - a - b1 * x[0][i] - b2 * x[1][i];

    }/*ELSE*/

  }/*THEN*/

  if (sd) {

    if (nrow == 0)
      *sd = R_NaN;
    else if (nrow <= 3)
      *sd = 0;
    else {

      if (resid)
        for (i = 0; i < nrow; i++)
          sse += resid[i] * resid[i];
      else if (fitted)
        for (i = 0; i < nrow; i++)
          sse += (y[i] - fitted[i]) * (y[i] - fitted[i]);
      else {

        if (singular1 && !singular2)
          for (i = 0; i < nrow; i++)
            sse += (y[i] - a - b2 * x[1][i]) * (y[i] - a - b2 * x[1][i]);
        else if (!singular1 && singular2)
          for (i = 0; i < nrow; i++)
            sse += (y[i] - a - b1 * x[0][i]) * (y[i] - a - b1 * x[0][i]);
        else if (singular1 && singular2)
          for (i = 0; i < nrow; i++)
            sse += (y[i] - a) * (y[i] - a);
        else
          for (i = 0; i < nrow; i++)
            sse += (y[i] - a - b1 * x[0][i] - b2 * x[1][i]) * (y[i] - a - b1 * x[0][i] - b2 * x[1][i]);

      }/*ELSE*/

      *sd = sqrt(sse / (nrow - 3));

    }/*ELSE*/

  }/*THEN*/

}/*C_OLS2*/

static void c_olsp(double **x, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd) {

double *qr = NULL;

  /* prepare the data for the QR decomposition. */
  qr = Calloc1D(nrow * (ncol + 1), sizeof(double));
  c_qr_matrix(qr, x, nrow, ncol, NULL, nrow);

  /* general case: multiple regression without missing data. */
  c_qr(qr, y, nrow, ncol + 1, fitted, resid, beta, sd);
  Free1D(qr);

}/*C_OLSP*/

/* fast implementation of least squares with just the intercept and with missing
 * data. */
static void c_ols0_with_missing(double *y, int nrow, double *fitted,
    double *resid, double *beta, double *sd) {

int i = 0, ncomplete = 0;
long double mean = 0, rsd = 0;

  /* scan the data to determine which observations are complete. */
  for (i = 0; i < nrow; i++) {

    if (!ISNAN(y[i])) {

      mean += y[i];
      ncomplete++;

    }/*THEN*/

  }/*FOR*/

  mean /= ncomplete;

  if (sd) {

    if (ncomplete == 0)
      *sd = R_NaN;
    else if (ncomplete == 1)
      *sd = 0;
    else {

      for (i = 0; i < nrow; i++)
        if (!ISNAN(y[i]))
          rsd += (y[i] - mean) * (y[i] - mean);

      *sd = sqrt(rsd / (ncomplete - 1));

    }/*ELSE*/

  }/*THEN*/

  /* the only coefficient is the intercept, which is the mean of the response. */
  if (beta)
    *beta = mean;
  if (fitted)
    for (i = 0; i < nrow; i++)
      fitted[i] = mean;
  if (resid)
    for (i  = 0; i < nrow; i++)
      resid[i] = y[i] - mean;

}/*C_OLS0_WITH_MISSING*/

/* general implementation of ordinary least squares with missing data. */
static void c_olsp_with_missing(double **x, double *y, int nrow, int ncol,
    double *fitted, double *resid, double *beta, double *sd) {

int i = 0, j = 0, k = 0, check = 0, ncomplete = 0;
int *complete = NULL;
double *qr = 0, *new_y = 0;

  /* scan the data to determine which observations are complete. */
  complete = Calloc1D(nrow, sizeof(int));

  for (i = 0; i < nrow; i++) {

    /* the response is missing. */
    check = !ISNAN(y[i]);

    /* one of the explanatory variables are missing. */
    for (j = 0; j < ncol; j++)
      check &= !ISNAN(x[j][i]);

    complete[i] = check;
    ncomplete += check;

  }/*FOR*/

  /* prepare the data for the QR decomposition. */
  qr = Calloc1D(ncomplete * (ncol + 1), sizeof(double));
  c_qr_matrix(qr, x, nrow, ncol, complete, ncomplete);
  new_y = Calloc1D(ncomplete, sizeof(double));
  for (i = 0, k = 0; i < nrow; i++)
    if (complete[i])
      new_y[k++] = y[i];
  /* general case: multiple regression with missing data. */
  c_qr(qr, new_y, ncomplete, ncol + 1, fitted, resid, beta, sd);

  /* rearrange the fitted values and the residuals to match the original data. */
  if (fitted) {

    for (i = nrow - 1, k = ncomplete - 1; i >= 0; i--) {

      if (complete[i])
        fitted[i] = fitted[k--];
      else
        fitted[i] = NA_REAL;

    }/*FOR*/

  }/*THEN*/

  if (resid) {

    for (i = nrow - 1, k = ncomplete - 1; i >= 0; i--) {

      if (complete[i])
        resid[i] = resid[k--];
      else
        resid[i] = NA_REAL;

    }/*FOR*/

  }/*THEN*/

  Free1D(complete);
  Free1D(new_y);
  Free1D(qr);

}/*C_OLSP_WITH_MISSING*/

/* compute least squares efficiently by special-casing whenever possible. */
void c_ols(double **x, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd, int missing) {

  if (!missing) {

    if (ncol == 0) {

      /* special case: null model. */
      c_ols0(y, nrow, fitted, resid, beta, sd);

    }/*THEN*/
    else if (ncol == 1) {

      /* special case: simple regression. */
      c_ols1(x, y, nrow, fitted, resid, beta, sd);

    }/*THEN*/
    else if (ncol == 2) {

      /* special case: two regression coefficients. */
      c_ols2(x, y, nrow, fitted, resid, beta, sd);

    }/*THEN*/
    else {

      /* general case: multiple regression without missing data. */
      c_olsp(x, y, nrow, ncol, fitted, resid, beta, sd);

    }/*ELSE*/

  }/*THEN*/
  else {

    if (ncol == 0) {

      /* special case: null model. */
      c_ols0_with_missing(y, nrow, fitted, resid, beta, sd);

    }/*THEN*/
    else {

      /* general case: multiple regression with missing data. */
      c_olsp_with_missing(x, y, nrow, ncol, fitted, resid, beta, sd);

    }/*ELSE*/

  }/*ELSE*/

}/*C_OLS*/

