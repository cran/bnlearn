#include "include/rcore.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/covariance.h"

/* -------------------- function declarations --------------------------- */

void c_svd(double *A, double *U, double *D, double *V, int *nrow,
    int *ncol, int *mindim, int strict, int *errcode);

/* -------------------- C level interfaces to LAPACK -------------------- */

/* C-level wrapper around the dgesdd() F77 routine. Note that the input
 * matrix A is overwritten by dgesdd(), so it's sensible to have a
 * backup copy in case it's needed later. */
void c_svd(double *A, double *U, double *D, double *V, int *nrow,
    int *ncol, int *mindim, int strict, int *errcode) {

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

/* helper function to allocate U, D and Vt for SVD. */
void c_udvt(double **u, double **d, double **vt, int ncol) {

  *u = Calloc1D(ncol * ncol, sizeof(double));
  *d = Calloc1D(ncol, sizeof(double));
  *vt = Calloc1D(ncol * ncol, sizeof(double));

}/*C_UDVT*/

/* C-level function to compute Moore-Penrose Generalized Inverse of a square matrix. */
void c_ginv(double *covariance, int ncol, double *mpinv) {

int i = 0, j = 0, errcode = 0;
double *u = NULL, *d = NULL, *vt = NULL, *backup = NULL;
double sv_tol = 0, zero = 0, one = 1;
char transa = 'N', transb = 'N';

  c_udvt(&u, &d, &vt, ncol);

  if (covariance != mpinv) {

    backup = Calloc1D(ncol * ncol, sizeof(double));
    memcpy(backup, covariance, ncol * ncol * sizeof(double));

  }/*THEN*/

  /* compute the SVD decomposition. */
  c_svd(covariance, u, d, vt, &ncol, &ncol, &ncol, FALSE, &errcode);

  /* if SVD fails, catch the error code and free all buffers. */
  if (errcode == 0) {

    /* set the threshold for the singular values as in corpcor. */
    sv_tol = ncol * d[0] * MACHINE_TOL * MACHINE_TOL;

    /* the first multiplication, U * D^{-1} is easy. */
    for (i = 0; i < ncol; i++)
      for (j = 0; j < ncol; j++)
        u[CMC(i, j, ncol)] = u[CMC(i, j, ncol)] * ((d[j] > sv_tol) ? 1/d[j] : 0);

    /* the second one, (U * D^{-1}) * Vt  is a real matrix multiplication. */
    F77_CALL(dgemm)(&transa, &transb, &ncol, &ncol, &ncol, &one, u,
      &ncol, vt, &ncol, &zero, mpinv, &ncol);

  }/*THEN*/

  if (covariance != mpinv) {

    memcpy(covariance, backup, ncol * ncol * sizeof(double));
    Free1D(backup);

  }/*THEN*/

  Free1D(u);
  Free1D(d);
  Free1D(vt);

  if (errcode)
    error("an error (%d) occurred in the call to c_ginv().\n", errcode);

}/*C_GINV*/

/* fast inverse for symmetric positive definite matries. */
void c_finv(double *cov, int *ncol, double *mpinv) {

double det = 0;
char u = 'U';
int i = 0, j = 0, errcode = 0;


  switch(*ncol) {

    case 2:

      det = cov[0] * cov[3] - cov[1] * cov[2];

      mpinv[0] =  cov[3] / det;
      mpinv[1] = -cov[2] / det;
      mpinv[2] = -cov[2] / det;
      mpinv[3] =  cov[0] / det;

      return;

    case 3:

      det = cov[0] * (cov[8] * cov[4] - cov[5] * cov[7]) -
            cov[1] * (cov[8] * cov[3] - cov[5] * cov[6]) +
            cov[2] * (cov[7] * cov[3] - cov[4] * cov[6]);

      mpinv[0] =  (cov[8] * cov[4] - cov[5] * cov[7]) / det;
      mpinv[1] = -(cov[8] * cov[1] - cov[2] * cov[7]) / det;
      mpinv[2] =  (cov[5] * cov[1] - cov[2] * cov[4]) / det;
      mpinv[3] = -(cov[8] * cov[3] - cov[5] * cov[6]) / det;
      mpinv[4] =  (cov[8] * cov[0] - cov[2] * cov[6]) / det;
      mpinv[5] = -(cov[5] * cov[0] - cov[2] * cov[1]) / det;
      mpinv[6] =  (cov[7] * cov[1] - cov[4] * cov[6]) / det;
      mpinv[7] = -(cov[7] * cov[0] - cov[1] * cov[6]) / det;
      mpinv[8] =  (cov[4] * cov[0] - cov[1] * cov[3]) / det;

      return;

    case 4:

      mpinv[0]  =  cov[5]  * cov[10] * cov[15] - cov[5]  * cov[11] * cov[14] -
                   cov[9]  * cov[6]  * cov[15] + cov[9]  * cov[7]  * cov[14] +
                   cov[13] * cov[6]  * cov[11] - cov[13] * cov[7]  * cov[10];
      mpinv[1]  = -cov[1]  * cov[10] * cov[15] + cov[1]  * cov[11] * cov[14] +
                   cov[9]  * cov[2]  * cov[15] - cov[9]  * cov[3]  * cov[14] -
                   cov[13] * cov[2]  * cov[11] + cov[13] * cov[3]  * cov[10];
      mpinv[2]  =  cov[1]  * cov[6]  * cov[15] - cov[1]  * cov[7]  * cov[14] -
                   cov[5]  * cov[2]  * cov[15] + cov[5]  * cov[3]  * cov[14] +
                   cov[13] * cov[2]  * cov[7]  - cov[13] * cov[3]  * cov[6];
      mpinv[3]  = -cov[1]  * cov[6]  * cov[11] + cov[1]  * cov[7]  * cov[10] +
                   cov[5]  * cov[2]  * cov[11] - cov[5]  * cov[3]  * cov[10] -
                   cov[9]  * cov[2]  * cov[7]  + cov[9]  * cov[3]  * cov[6];
      mpinv[4]  = -cov[4]  * cov[10] * cov[15] + cov[4]  * cov[11] * cov[14] +
                   cov[8]  * cov[6]  * cov[15] - cov[8]  * cov[7]  * cov[14] -
                   cov[12] * cov[6]  * cov[11] + cov[12] * cov[7]  * cov[10];
      mpinv[5]  =  cov[0]  * cov[10] * cov[15] - cov[0]  * cov[11] * cov[14] -
                   cov[8]  * cov[2]  * cov[15] + cov[8]  * cov[3]  * cov[14] +
                   cov[12] * cov[2]  * cov[11] - cov[12] * cov[3]  * cov[10];
      mpinv[6]  = -cov[0]  * cov[6]  * cov[15] + cov[0]  * cov[7]  * cov[14] +
                   cov[4]  * cov[2]  * cov[15] - cov[4]  * cov[3]  * cov[14] -
                   cov[12] * cov[2]  * cov[7]  + cov[12] * cov[3]  * cov[6];
      mpinv[7]  =  cov[0]  * cov[6]  * cov[11] - cov[0]  * cov[7]  * cov[10] -
                   cov[4]  * cov[2]  * cov[11] + cov[4]  * cov[3]  * cov[10] +
                   cov[8]  * cov[2]  * cov[7]  - cov[8]  * cov[3]  * cov[6];
      mpinv[8]  =  cov[4]  * cov[9]  * cov[15] - cov[4]  * cov[11] * cov[13] -
                   cov[8]  * cov[5]  * cov[15] + cov[8]  * cov[7]  * cov[13] +
                   cov[12] * cov[5]  * cov[11] - cov[12] * cov[7]  * cov[9];
      mpinv[9]  = -cov[0]  * cov[9]  * cov[15] + cov[0]  * cov[11] * cov[13] +
                   cov[8]  * cov[1]  * cov[15] - cov[8]  * cov[3]  * cov[13] -
                   cov[12] * cov[1]  * cov[11] + cov[12] * cov[3]  * cov[9];
      mpinv[11] = -cov[0]  * cov[5]  * cov[11] + cov[0]  * cov[7]  * cov[9] +
                   cov[4]  * cov[1]  * cov[11] - cov[4]  * cov[3]  * cov[9] -
                   cov[8]  * cov[1]  * cov[7]  + cov[8]  * cov[3]  * cov[5];
      mpinv[10] =  cov[0]  * cov[5]  * cov[15] - cov[0]  * cov[7]  * cov[13] -
                   cov[4]  * cov[1]  * cov[15] + cov[4]  * cov[3]  * cov[13] +
                   cov[12] * cov[1]  * cov[7]  - cov[12] * cov[3]  * cov[5];
      mpinv[12] = -cov[4]  * cov[9]  * cov[14] + cov[4]  * cov[10] * cov[13] +
                   cov[8]  * cov[5]  * cov[14] - cov[8]  * cov[6]  * cov[13] -
                   cov[12] * cov[5]  * cov[10] + cov[12] * cov[6]  * cov[9];
      mpinv[13] =  cov[0]  * cov[9]  * cov[14] - cov[0]  * cov[10] * cov[13] -
                   cov[8]  * cov[1]  * cov[14] + cov[8]  * cov[2]  * cov[13] +
                   cov[12] * cov[1]  * cov[10] - cov[12] * cov[2]  * cov[9];
      mpinv[14] = -cov[0]  * cov[5]  * cov[14] + cov[0]  * cov[6]  * cov[13] +
                   cov[4]  * cov[1]  * cov[14] - cov[4]  * cov[2]  * cov[13] -
                   cov[12] * cov[1]  * cov[6]  + cov[12] * cov[2]  * cov[5];
      mpinv[15] =  cov[0]  * cov[5]  * cov[10] - cov[0]  * cov[6]  * cov[9] -
                   cov[4]  * cov[1]  * cov[10] + cov[4]  * cov[2]  * cov[9] +
                   cov[8]  * cov[1]  * cov[6]  - cov[8]  * cov[2]  * cov[5];

      det = cov[0] * mpinv[0] + cov[1] * mpinv[4] + cov[2] * mpinv[8] + cov[3] * mpinv[12];

      for (i = 0; i < 16; i++)
        mpinv[i] = mpinv[i] / det;

      return;

    default:
      break;

  }/*SWITCH*/

  /* copy the original matrix, it gets overwritten otherwise. */
  memcpy(mpinv, cov, (*ncol) * (*ncol) * sizeof(double));

  /* compute the upper triangular part of the inverse. */
  F77_CALL(dpotrf)(&u, ncol, mpinv, ncol, &errcode);
  F77_CALL(dpotri)(&u, ncol, mpinv, ncol, &errcode);

  /* fill in the lower trinagular part of the matrix. */
  for (i = 0; i < *ncol; i++)
    for (j = i + 1; j < *ncol; j++)
      mpinv[CMC(j, i, *ncol)] = mpinv[CMC(i, j, *ncol)];

}/*C_FINV*/

/* C-level function to compute the quadratic form (xT)S^{-1}(y). */
double c_quadratic(double *x, int *ncol, double *sigma, double *y, double *workspace) {

char transa = 'N', transb = 'N';
double zero = 0, d_one = 1, res = 0;
int i = 0, i_one = 1;

  /* first compute the (xT)S^{-1} part ... */
  F77_CALL(dgemm)(&transa, &transb, &i_one, ncol, ncol, &d_one, x,
    &i_one, sigma, ncol, &zero, workspace, &i_one);

  /* ... then multiply it again for x. */
  for (i = 0; i < *ncol; i++)
    res += workspace[i] * y[i];

  return res;

}/*C_QUADRATIC*/

/* C-level function to compute S1 * (S2 * x + a * mu). */
void c_rotate(double *S1, double *S2, double *x, double *a, double *mu,
    int *ncol, double *workspace) {

char transa = 'N', transb = 'N';
double zero = 0, d_one = 1;
int i_one = 1;

  /* this double call to dgemm() requires a save-restore for mu. */
  memcpy(workspace, mu, *ncol * sizeof(double));

  F77_CALL(dgemm)(&transa, &transb, ncol, &i_one, ncol, &d_one, S2,
    ncol, x, ncol, a, mu, ncol);

  F77_CALL(dgemm)(&transa, &transb, ncol, &i_one, ncol, &d_one, S1,
    ncol, mu, ncol, &zero, x, ncol);

  memcpy(mu, workspace, *ncol * sizeof(double));

}/*C_ROTATE*/

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

  if (info != 0)
    error("an error (%d) occurred in the call to dqrsl().\n", &info);

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
  c_sd(rsd, nrow, ncol, 0, FALSE, sd);

  if (!resid)
    Free1D(rsd);
  if (!fitted)
    Free1D(ftt);
  Free1D(pivot);
  Free1D(work);
  Free1D(qraux);

}/*C_QR*/

