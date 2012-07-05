#include "common.h"
#include <Rmath.h>
#include <R_ext/Lapack.h>

/* -------------------- function declarations --------------------------- */

void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim, int strict, int *errcode);
double c_det(double *matrix, int *rows);

/* -------------------- R level interfaces to LAPACK -------------------- */

/* Determinant. */
SEXP r_det(SEXP matrix, int scale) {

int i = 0, nr = nrows(matrix);
short int duplicated = 0;
SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));

  /* duplicate the matrix so the original one is not overwritten. */
  if ((duplicated = NAMED(matrix)) > 0)
    PROTECT(matrix = duplicate(matrix));

  /* rescale the determinant by multiplying the matrix. */
  for (i = 0; i < nr * nr; i++)
    REAL(matrix)[i] *= scale;

  /* call Lapack to do the grunt work. */
  NUM(result) = c_det(REAL(matrix), &nr);

  if (duplicated > 0)
    UNPROTECT(2);
  else
    UNPROTECT(1);

  return result;

}/*R_DET*/

/* -------------------- C level interfaces to LAPACK -------------------- */

/* C-level wrapper around the dgesvd() F77 routine. Note that the input
 * matrix A is overwritten by dgesvd(), so it's sensible to have a
 * backup copy in case it's needed later. */
void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim, int strict, int *errcode) {

int lwork = -1, *iwork = NULL;
char jobz = 'A';
double tmp = 0, *work = NULL;

  iwork = (int *)Calloc(8 * (*mindim), int);

  /* ask for the optimal size of the work array. */
  F77_CALL(dgesdd)(&jobz, nrows, ncols, A, nrows, D, U, nrows,
                   V, mindim, &tmp, &lwork, iwork, errcode);

  lwork = (int)tmp;
  work = (double *)Calloc(lwork, double);

  /* actual call */
  F77_NAME(dgesdd)(&jobz, nrows, ncols, A, nrows, D, U, nrows,
                   V, mindim, work, &lwork, iwork, errcode);

  Free(work);
  Free(iwork);

  if (*errcode && strict)
    error("an error (%d) occurred in the call to dgesdd().\n", *errcode);

}/*C_SVD*/

/* C-level function to compute the determinant of a real-valued square matrix,
 * modeled after the moddet_ge_real() function in Lapack.c. */
double c_det(double *matrix, int *rows) {

int sign = 1, i = 0, info = 0, *jpvt = NULL;
double det = 1;

  jpvt = (int *) Calloc(*rows, int);

  /* comute the A = L*U decomposition. */
  F77_CALL(dgetrf)(rows, rows, matrix, rows, jpvt, &info);

  if (info < 0) {

    error("an error (%d) occurred in the call to dgesvd().\n", info);

  }/*THEN*/
  else if (info > 0) {

    /* the matrix is singular, so the determinant is zero. */
    det = 0;

  }/*THEN*/
  else {

    /* the matrix is full rank, compute the determinant. */
    for (i = 0; i < *rows; i++) {

      if (jpvt[i] != (i + 1))
        sign = -sign;

      det *= matrix[CMC(i, i, *rows)];

    }/*FOR*/

  }/*ELSE*/

  Free(jpvt);

  return sign * det;

}/*C_DET*/

/* C-level function to compute Moore-Penrose Generalized Inverse of a square matrix. */
void c_ginv(double *covariance, int *ncols, double *mpinv) {

int i = 0, j = 0, errcode = 0;
double *u = NULL, *d = NULL, *vt = NULL, *backup = NULL;
double tol = MACHINE_TOL, zero = 0, one = 1;
char transa = 'N', transb = 'N';

  u = Calloc((*ncols) * (*ncols), double);
  memset(u, '\0', (*ncols) * (*ncols) * sizeof(double));
  d = Calloc((*ncols), double);
  memset(d, '\0', (*ncols) * sizeof(double));
  vt = Calloc((*ncols) * (*ncols), double);
  memset(vt, '\0', (*ncols) * (*ncols) * sizeof(double));

  if (covariance != mpinv) {

    backup = Calloc((*ncols) * (*ncols), double);
    memcpy(backup, covariance, (*ncols) * (*ncols) * sizeof(double));

  }/*THEN*/

  /* compute the SVD decomposition. */
  c_svd(covariance, u, d, vt, ncols, ncols, ncols, FALSE, &errcode);

  /* if SVD fails, catch the error code and free all buffers. */
  if (errcode)
    goto end;

  /* the first multiplication, U * D^{-1} is easy. */
  for (i = 0; i < *ncols; i++)
    for (j = 0; j < *ncols; j++)
      u[CMC(i, j, *ncols)] = u[CMC(i, j, *ncols)] * ((d[j] > tol) ? 1/d[j] : 0);

  /* the second one, (U * D^{-1}) * Vt  is a real matrix multiplication. */
  F77_CALL(dgemm)(&transa, &transb, ncols, ncols, ncols, &one, u,
    ncols, vt, ncols, &zero, mpinv, ncols);

end:

  if (covariance != mpinv) {

    memcpy(covariance, backup, (*ncols) * (*ncols) * sizeof(double));
    Free(backup);

  }/*THEN*/

  Free(u);
  Free(d);
  Free(vt);

  if (errcode)
    error("an error (%d) occurred in the call to c_ginv().\n", errcode);

}/*C_GINV*/

/* fast inverse for symmetric positive definite matries. */
void c_finv(double *cov, int *ncols, double *mpinv) {

double det = 0;
char u = 'U';
int i = 0, j = 0, errcode = 0;


  switch(*ncols) {

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
  memcpy(mpinv, cov, (*ncols) * (*ncols) * sizeof(double));

  /* compute the upper triangular part of the inverse. */
  F77_CALL(dpotrf)(&u, ncols, mpinv, ncols, &errcode);
  F77_CALL(dpotri)(&u, ncols, mpinv, ncols, &errcode);

  /* fill in the lower trinagular part of the matrix. */
  for (i = 0; i < *ncols; i++)
    for (j = i + 1; j < *ncols; j++)
      mpinv[CMC(j, i, *ncols)] = mpinv[CMC(i, j, *ncols)];

}/*C_FINV*/

/* C-level function to compute the quadratic form (xT)S^{-1}(y). */
double c_quadratic(double *x, int *ncols, double *sigma, double *y, double *workspace) {

char transa = 'N', transb = 'N';
double zero = 0, d_one = 1, res = 0;
int i = 0, i_one = 1;

  /* first compute the (xT)S^{-1} part ... */
  F77_CALL(dgemm)(&transa, &transb, &i_one, ncols, ncols, &d_one, x,
    &i_one, sigma, ncols, &zero, workspace, &i_one);

  /* ... then multiply it again for x. */
  for (i = 0; i < *ncols; i++)
    res += workspace[i] * y[i];

  return res;

}/*C_QUADRATIC*/

/* C-level function to compute S1 * (S2 * x + a * mu). */
void c_rotate(double *S1, double *S2, double *x, double *a, double *mu,
    int *ncols, double *workspace) {

char transa = 'N', transb = 'N';
double zero = 0, d_one = 1;
int i_one = 1;

  /* this double call to dgemm() requires a save-restore for mu. */
  memcpy(workspace, mu, *ncols * sizeof(double));

  F77_CALL(dgemm)(&transa, &transb, ncols, &i_one, ncols, &d_one, S2,
    ncols, x, ncols, a, mu, ncols);

  F77_CALL(dgemm)(&transa, &transb, ncols, &i_one, ncols, &d_one, S1,
    ncols, mu, ncols, &zero, x, ncols);

  memcpy(mu, workspace, *ncols * sizeof(double));

}/*C_ROTATE*/

