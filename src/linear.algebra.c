#include "include/rcore.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/covariance.h"

/* -------------------- C level interfaces to LAPACK -------------------- */

/* C-level function to compute the log-determinant of a symmetric real matrix,
   which can be either positive or NaN. */
double c_logdet(double *matrix, int rows) {

int sign = 1, i = 0, info = 0, *jpvt = NULL;
long double logdet = 0;

  jpvt = (int *) Calloc(rows, int);

  /* compute the A = L*U decomposition. */
  F77_CALL(dgetrf)(&rows, &rows, matrix, &rows, jpvt, &info);

  if (info < 0) {

    error("an error (%d) occurred in the call to dgesvd().\n", info);

  }/*THEN*/
  else if (info > 0) {

    /* the matrix is singular, so the determinant is zero. */
    logdet = R_NegInf;

  }/*THEN*/
  else {

    /* the matrix is full rank, compute the determinant. */
    for (i = 0; i < rows; i++) {

      /* cumulate the sign changes introduced by pivoting... */
      if (jpvt[i] != (i + 1))
        sign = -sign;

      /* ... and the signs of the eigenvalues. */
      if (matrix[CMC(i, i, rows)] < 0)
        sign = -sign;

      /* compute the determinant on the log-scale disregarding the sign. */
      logdet += log(fabs(matrix[CMC(i, i, rows)]));

    }/*FOR*/

  }/*ELSE*/

  Free(jpvt);

  return (sign > 0) ? (double) logdet : R_NaN;

}/*C_LOGDET*/

/* C-level wrapper around the dgesdd() F77 routine. Note that the input
 * matrix A is overwritten by dgesdd(), so it's sensible to have a
 * backup copy in case it's needed later. */
void c_svd(double *A, double *U, double *D, double *V, int *nrow,
    int *ncol, int *mindim, bool strict, int *errcode) {

int lwork = -1, *iwork = NULL;
char jobz = 'A';
double tmp = 0, *work = NULL;

  iwork = Calloc1D(8 * (*mindim), sizeof(int));

  /* ask for the optimal size of the work array. */
  F77_CALL(dgesdd)(&jobz, nrow, ncol, A, nrow, D, U, nrow,
                   V, mindim, &tmp, &lwork, iwork, errcode);

  lwork = (int)tmp;
  work = Calloc1D(lwork, sizeof(double));

  /* actual call */
  F77_NAME(dgesdd)(&jobz, nrow, ncol, A, nrow, D, U, nrow,
                   V, mindim, work, &lwork, iwork, errcode);

  Free1D(work);
  Free1D(iwork);

  if (*errcode && strict)
    error("an error (%d) occurred in the call to dgesdd().\n", *errcode);

}/*C_SVD*/

/* C-level function to compute Moore-Penrose Generalized Inverse of a square matrix. */
void c_ginv(covariance cov, covariance mpinv) {

int i = 0, j = 0, errcode = 0;
double sv_tol = 0, zero = 0, one = 1;
char transa = 'N', transb = 'N';
covariance backup = (covariance){ 0 };

  if (cov.mat != mpinv.mat) {

    backup = new_covariance(cov.dim, TRUE);
    copy_covariance(&cov, &backup);

  }/*THEN*/

  /* compute the SVD decomposition. */
  c_svd(cov.mat, cov.u, cov.d, cov.vt, &cov.dim, &cov.dim, &cov.dim,
    FALSE, &errcode);

  /* if SVD fails, catch the error code and free all buffers. */
  if (errcode == 0) {

    /* set the threshold for the singular values as in corpcor. */
    sv_tol = cov.dim * cov.d[0] * MACHINE_TOL * MACHINE_TOL;

    /* the first multiplication, U * D^{-1} is easy. */
    for (i = 0; i < cov.dim; i++)
      for (j = 0; j < cov.dim; j++)
        cov.u[CMC(i, j, cov.dim)] = cov.u[CMC(i, j, cov.dim)] *
                                      ((cov.d[j] > sv_tol) ? 1/cov.d[j] : 0);

    /* the second one, (U * D^{-1}) * Vt  is a real matrix multiplication. */
    F77_CALL(dgemm)(&transa, &transb, &cov.dim, &cov.dim, &cov.dim, &one, cov.u,
      &cov.dim, cov.vt, &cov.dim, &zero, mpinv.mat, &cov.dim);

  }/*THEN*/

  if (cov.mat != mpinv.mat) {

    copy_covariance(&backup, &cov);
    FreeCOV(backup);

  }/*THEN*/

  if (errcode)
    error("an error (%d) occurred in the call to c_ginv().\n", errcode);

}/*C_GINV*/

/* C-level function to perform OLS via QR decomposition. */
void c_qr(double *qr, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd) {

int i = 0, job = 10, rank = 0, info = 0, *pivot = NULL, pivoted = FALSE;
double tol = MACHINE_TOL, *qraux = NULL, *work = NULL;
double *bb = NULL, *rsd = NULL, *ftt = NULL;

  /* special case for sample size = 1. */
  if (nrow == 1) {

    if (beta) {

      beta[0] = *y;
      for (i = 1; i < ncol; i++)
        beta[i] = NA_REAL;

    }/*FIT*/

    if (fitted)
      *fitted = *y;

    if (resid)
      *resid = 0;

    *sd = 0;

    return;

  }/*THEN*/

  /* allocate the working space. */
  qraux = Calloc1D(ncol, sizeof(double));
  work = Calloc1D(2 * ncol, sizeof(double));
  pivot = Calloc1D(ncol, sizeof(int));
  for (i = 0; i < ncol; i++)
    pivot[i] = i + 1;

  /* set the job magic number depending on the targets of the computation. */
  if (fitted) {

    job += 1;
    ftt = fitted;

  }/*THEN*/
  else {

    ftt = Calloc1D(nrow, sizeof(double));

  }/*ELSE*/
  if (beta) {

    job += 100;
    bb = Calloc1D(ncol, sizeof(double));

  }/*THEN*/
  if (resid)
    rsd = resid;
  else
    rsd = Calloc1D(nrow, sizeof(double));

  /* perform the QR decomposition. */
  F77_CALL(dqrdc2)(qr, &nrow, &nrow, &ncol, &tol, &rank, qraux, pivot, work);

  /* operate on a backup copy of the response variable. */
  memcpy(ftt, y, nrow * sizeof(double));

  /* compute the fitted values, residuals, coefficients and standard errors. */
  /*       dqrsl( x,  ldx,   n,     k,     qraux, y,   qy,     qty, */
  F77_CALL(dqrsl)(qr, &nrow, &nrow, &rank, qraux, ftt, NULL,   ftt,
  /*  b,      rsd,    xb,  job,  info) */
      bb,   rsd,    ftt, &job, &info);

  /* 'info' is the error code, if it is different from zero something went
   * wrong and the model is perfectly singular to the point the QR
   * decomposition cannot be computed at all. */
  if (info != 0) {

    /* set all coefficients to NA. */
    if (beta)
      for (i = 0; i < ncol; i++)
        beta[i] = NA_REAL;

    /* set all residuals to zero. */
    memset(rsd, '\0', nrow * sizeof(double));
    /* set the standard error to zero. */
    *sd = 0;
    /* set all the fitted values to the corresponding values in the response. */
    if (fitted)
      memcpy(ftt, y, nrow * sizeof(double));

  }/*THEN*/
  else {

    if (beta) {

      /* set the coefficients to NA in rank-deficient problems; they are moved
       * to the end of the array. */
      if (rank < ncol)
        for (i = rank; i < ncol; i++)
          bb[i] = NA_REAL;

      /* check whether the coefficients have been pivoted (R does that as a
       * separate check from the rank check, so let's do the same). */
      for (i = 0; i < ncol; i++)
        if (pivot[i] != i + 1) {

          pivoted = TRUE;
          break;

        }/*THEN*/

      /* coefficients are pivoted in singular problems, move them back. */
      if (pivoted) {

        for (i = 0; i < ncol; i++)
          beta[pivot[i] - 1] = bb[i];

      }/*THEN*/
      else {

        memcpy(beta, bb, ncol * sizeof(double));

      }/*ELSE*/

      Free1D(bb);

    }/*THEN*/

    /* compute the standard deviation of the residuals. */
    c_sd(rsd, nrow, ncol, 0, sd);

  }/*ELSE*/

  if (!resid)
    Free1D(rsd);
  if (!fitted)
    Free1D(ftt);
  Free1D(pivot);
  Free1D(work);
  Free1D(qraux);

}/*C_QR*/

