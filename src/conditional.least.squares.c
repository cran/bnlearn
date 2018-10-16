#include "include/rcore.h"
#include "include/covariance.h"
#include "include/blas.h"
#include "include/globals.h"
#include "include/matrix.h"

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

      cc[i * 2 + 1] = (fabsl(var_x[i]) < MACHINE_TOL) ? 0 : cov[i] / var_x[i];
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

    for (i = 0; i < ncond; i++) {

      if (nz[i] == 0)
        sd[i] = R_NaN;
      else if (nz[i] <= 2)
        sd[i] = 0;
      else
        sd[i] = sqrt(ssr[i] / (nz[i] - 2));

    }/*FOR*/

    Free1D(ssr);

  }/*THEN*/

  Free1D(cc);
  Free1D(mean_y);
  Free1D(mean_x);
  Free1D(var_x);
  Free1D(cov);
  Free1D(nz);

}/*C_CLS1*/

/* fast implementation of conditional least squares with two regression
 * coefficients. */
static void c_cls2(double **x, double *y, int *z, int nrow, int ncond,
    double *fitted, double *resid, double *beta, double *sd) {

int i = 0, *nz = NULL, singular1 = FALSE, singular2 = FALSE;
double *cc = NULL;
long double *mean_y = NULL, *mean_x1 = NULL, *mean_x2 = NULL;
long double **var = NULL, **cov = NULL;
long double den = 0, *ssr = NULL;

  /* first pass: compute the conditional means of y and x1. */
  mean_y = Calloc1D(ncond, sizeof(long double));
  mean_x1 = Calloc1D(ncond, sizeof(long double));
  mean_x2 = Calloc1D(ncond, sizeof(long double));
  var = (long double **) Calloc2D(ncond, 3, sizeof(long double));
  cov = (long double **) Calloc2D(ncond, 3, sizeof(long double));
  nz = Calloc1D(ncond, sizeof(int));
  cc = Calloc1D(3 * ncond, sizeof(double));

  for (i = 0; i < nrow; i++) {

    mean_y[z[i] - 1] += y[i];
    mean_x1[z[i] - 1] += x[0][i];
    mean_x2[z[i] - 1] += x[1][i];
    nz[z[i] - 1]++;

  }/*FOR*/

  for (i = 0; i < ncond; i++) {

    mean_y[i] /= nz[i];
    mean_x1[i] /= nz[i];
    mean_x2[i] /= nz[i];

  }/*FOR*/

  /* second pass: compute the variances of (y, x1, x2); and the covariances of
   * ((y, x1), (y, x2), (x1, x2)). */
  for (i = 0; i < nrow; i++) {

    var[z[i] - 1][0] += (y[i] - mean_y[z[i] - 1]) * (y[i] - mean_y[z[i] - 1]);
    var[z[i] - 1][1] += (x[0][i] - mean_x1[z[i] - 1]) * (x[0][i] - mean_x1[z[i] - 1]);
    var[z[i] - 1][2] += (x[1][i] - mean_x2[z[i] - 1]) * (x[1][i] - mean_x2[z[i] - 1]);
    cov[z[i] - 1][0] += (y[i] - mean_y[z[i] - 1]) * (x[0][i] - mean_x1[z[i] - 1]);
    cov[z[i] - 1][1] += (y[i] - mean_y[z[i] - 1]) * (x[1][i] - mean_x2[z[i] - 1]);
    cov[z[i] - 1][2] += (x[0][i] - mean_x1[z[i] - 1]) * (x[1][i] - mean_x2[z[i] - 1]);

  }/*FOR*/

  /* the regression coefficients are computed using the closed form estimates
   * for the two-variable regression, which are effectively the same as the
   * corresponding partial correlation estimates. */
  for (i = 0; i < ncond; i++) {

    if (nz[i] == 0) {

      cc[i * 3] = cc[i * 3 + 1] = cc[i * 3 + 2] = R_NaN;

    }/*THEN*/
    else {

      /* there are three possible collinear configurations:
       *   1) the first variable is constant and collinear with the response;
       *   2) the second variable is constant and collinear with the response;
       *   3) the two variables are collinear with each other. */
      singular1 = (fabsl(var[i][1]) < MACHINE_TOL);
      singular2 = (fabsl(var[i][2]) < MACHINE_TOL) ||
                  (fabsl(cov[i][2]) / sqrt(var[i][1] * var[i][2]) > 1 - MACHINE_TOL);

      if (singular1 && !singular2) {

        cc[i * 3 + 1] = 0;
        cc[i * 3 + 2] = cov[i][1] / var[i][2];
        cc[i * 3] = mean_y[i] - mean_x2[i] * cc[i * 3 + 2];

      }/*THEN*/
      else if (!singular1 && singular2) {

        cc[i * 3 + 1] = cov[i][0] / var[i][1];
        cc[i * 3 + 2] = 0;
        cc[i * 3] = mean_y[i] - mean_x1[i] * cc[i * 3 + 1];

      }/*THEN*/
      else if (singular1 && singular2) {

       cc[i * 3 + 1] = cc[i * 3 + 2] = 0;
       cc[i * 3] = mean_y[i];

      }/*THEN*/
      else {

        den = (var[i][1] * var[i][2] - cov[i][2] * cov[i][2]);

        cc[i * 3 + 1] = (var[i][2] * cov[i][0] - cov[i][2] * cov[i][1]) / den;
        cc[i * 3 + 2] = (var[i][1] * cov[i][1] - cov[i][2] * cov[i][0]) / den;
        cc[i * 3] = mean_y[i] - mean_x1[i] * cc[i * 3 + 1] - mean_x2[i] * cc[i * 3 + 2];

      }/*ELSE*/

    }/*ELSE*/

  }/*FOR*/

  if (beta)
    memcpy(beta, cc, 3 * ncond * sizeof(double));
  if (fitted)
    for (i  = 0; i < nrow; i++)
      fitted[i] = cc[(z[i] - 1) * 3] + x[0][i] * cc[(z[i] - 1) * 3 + 1]
                    + x[1][i] * cc[(z[i] - 1) * 3 + 2];
  if (resid) {

    if (fitted)
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - fitted[i];
    else
      for (i  = 0; i < nrow; i++)
        resid[i] = y[i] - cc[(z[i] - 1) * 3] - x[0][i] * cc[(z[i] - 1) * 3 + 1]
                     - x[1][i] * cc[(z[i] - 1) * 3 + 2];

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
        ssr[z[i] - 1] +=
          (y[i] - cc[(z[i] - 1) * 3] - x[0][i] * cc[(z[i] - 1) * 3 + 1]
                - x[1][i] * cc[(z[i] - 1) * 3 + 2]) *
          (y[i] - cc[(z[i] - 1) * 3] - x[0][i] * cc[(z[i] - 1) * 3 + 1]
                - x[1][i] * cc[(z[i] - 1) * 3 + 2]);

    for (i = 0; i < ncond; i++) {

      if (nz[i] == 0)
        sd[i] = R_NaN;
      else if (nz[i] <= 3)
        sd[i] = 0;
      else
        sd[i] = sqrt(ssr[i] / (nz[i] - 3));

    }/*FOR*/

    Free1D(ssr);

  }/*THEN*/

  Free1D(mean_y);
  Free1D(mean_x1);
  Free1D(mean_x2);
  Free2D(var, ncond);
  Free2D(cov, ncond);
  Free1D(nz);
  Free1D(cc);

}/*C_CLS2*/

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
      if (ISNAN(beta[i]))
        beta[i] = 0;

  Free1D(nz);
  Free1D(counter);
  Free2D(yy, ncond);
  Free2D(qr_matrix, ncond);
  if (fitted || resid)
    Free2D(zid, ncond);

}/*C_CLSP*/

/* fast implementation of conditional least squares with just the intercept. */
static void c_cls0_with_missing(double *y, int *z, int nrow, int ncond,
    double *fitted, double *resid, double *beta, double *sd) {

int i = 0, *nz = NULL;
long double *means = NULL, *ssr = NULL;

  /* the matrix of regression coefficients contains the conditional mean given
   * each configuration; estimate in a single pass (set to NaN if there a no
   * observations for a particular observation). */
  means = Calloc1D(ncond, sizeof(long double));
  nz = Calloc1D(ncond, sizeof(int));

  for (i = 0; i < nrow; i++) {

    if (!ISNAN(y[i]) && (z[i] != NA_INTEGER)) {

      means[z[i] - 1] += y[i];
      nz[z[i] - 1]++;

    }/*THEN*/

  }/*FOR*/

  for (i = 0; i < ncond; i++) {

    if (nz[i] == 0)
      means[i] = R_NaN;
    else
      means[i] /= nz[i];

  }/*FOR*/

  if (beta)
    for (i = 0; i < ncond; i++)
      beta[i] = means[i];
  if (fitted)
    for (i = 0; i < nrow; i++)
      fitted[i] = (z[i] != NA_INTEGER) && !ISNAN(y[i]) ? means[z[i] - 1] : NA_REAL;
  if (resid)
    for (i = 0; i < nrow; i++)
      resid[i] = (z[i] != NA_INTEGER) && !ISNAN(y[i]) ? y[i] - means[z[i] - 1] : NA_REAL;
  if (sd) {

    ssr = Calloc1D(ncond, sizeof(long double));

    if (resid) {

      for (i = 0; i < nrow; i++)
        if (!ISNAN(resid[i]))
          ssr[z[i] - 1] += resid[i] * resid[i];

    }/*THEN*/
    else if (fitted) {

      for (i = 0; i < nrow; i++)
        if (!ISNAN(fitted[i]))
          ssr[z[i] - 1] += (y[i] - fitted[i]) * (y[i] - fitted[i]);

    }/*THEN*/
    else {

      for (i = 0; i < nrow; i++)
        if ((z[i] != NA_INTEGER) && !ISNAN(y[i]))
          ssr[z[i] - 1] += (y[i] - means[z[i] - 1]) * (y[i] - means[z[i] - 1]);

    }/*ELSE*/

    for (i = 0; i < ncond; i++) {

      if (nz[i] == 0)
        sd[i] = R_NaN;
      else if (nz[i] == 1)
        sd[i] = 0;
      else
        sd[i] = sqrt(ssr[i] / (nz[i] - 1));

    }/*FOR*/

    Free1D(ssr);

  }/*THEN*/

  Free1D(nz);
  Free1D(means);

}/*C_CLS0_WITH_MISSING*/

/* general implementation of conditional least squares. */
void c_clsp_with_missing(double **x, double *y, int *z, int nrow, int ncol,
    int ncond, double *fitted, double *resid, double *beta, double *sd) {

int i = 0, j = 0, *nz = NULL, *counter = NULL, qr_ncol = ncol + 1;
double **qr_matrix = NULL, **yy = NULL, *ff = NULL, *rr = NULL;
int check = FALSE, ncomplete = 0, *complete = NULL, **zid = NULL;

  /* first pass: count how many observations for each conditioning. */
  nz = Calloc1D(ncond, sizeof(int));
  complete = Calloc1D(nrow, sizeof(int));

  for (i = 0; i < nrow; i++) {

    /* the response or the conditioning variables are missing. */
    check = !ISNAN(y[i]) && (z[i] != NA_INTEGER);

    /* one of the explanatory variables are missing. */
    for (j = 0; j < ncol; j++)
      check &= !ISNAN(x[j][i]);

    complete[i] = check;
    ncomplete += check;
    if (check)
      nz[z[i] - 1]++;

  }/*FOR*/

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

    if (!complete[i])
      continue;

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
      if (complete[i])
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

  /* make sure to put missing values into fitted values and residuals in the
   * positions corresponding to incomplete observations.*/
  if (fitted) {

    for (i = 0; i < nrow; i++)
      if (!complete[i])
        fitted[i] = NA_REAL;

  }/*THEN*/

  if (resid) {

    for (i = 0; i < nrow; i++)
      if (!complete[i])
        resid[i] = NA_REAL;

  }/*THEN*/

  /* post-process the regression coefficients and replace NAs with zeros, while
   * leaving NaNs in place for unobserved conditionings. */
  if (beta)
    for (i = 0; i < ncond * qr_ncol; i++)
      if (ISNAN(beta[i]))
        beta[i] = 0;

  Free1D(nz);
  Free1D(counter);
  Free1D(complete);
  Free2D(yy, ncond);
  Free2D(qr_matrix, ncond);
  if (fitted || resid)
    Free2D(zid, ncond);

}/*C_CLSP_WITH_MISSING*/

/* compute conditional least squares efficiently by special-casing whenever
 * possible. */
void c_cls(double **x, double *y, int *z, int nrow, int ncol, int ncond,
    double *fitted, double *resid, double *beta, double *sd, int missing) {

  if (!missing) {

    if (ncol == 0) {

      c_cls0(y, z, nrow, ncond, fitted, resid, beta, sd);

    }/*THEN*/
    else if (ncol == 1) {

      c_cls1(x, y, z, nrow, ncond, fitted, resid, beta, sd);

    }/*THEN*/
    else if (ncol == 2) {

      c_cls2(x, y, z, nrow, ncond, fitted, resid, beta, sd);

    }/*THEN*/
    else {

      c_clsp(x, y, z, nrow, ncol, ncond, fitted, resid, beta, sd);

    }/*ELSE*/

  }/*THEN*/
  else {

    if (ncol == 0) {

      c_cls0_with_missing(y, z, nrow, ncond, fitted, resid, beta, sd);

    }/*THEN*/
    else {

      c_clsp_with_missing(x, y, z, nrow, ncol, ncond, fitted, resid, beta, sd);

    }/*ELSE*/

  }/*ELSE*/

}/*C_CLS*/


