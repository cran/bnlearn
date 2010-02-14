#include "common.h"
#include <Rmath.h>
#include <R_ext/Lapack.h>

/* -------------------- function declarations --------------------------- */

static void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim);
double c_det(double *matrix, int *rows);

/* -------------------- R level interfaces to LAPACK -------------------- */

/* Singular Value Decomposition (SVD). */
SEXP r_svd(SEXP matrix) {

SEXP U, D, Vt, result, elnames;
int nr = 0, nc = 0, mindim = 0;
short int duplicated = 0;

  /* compute the dimensions of the SVD components. */
  nr = nrows(matrix);
  nc = ncols(matrix);
  mindim = imin2(nr, nc);

  /* allocate the U, D and t(V) matrices. */
  PROTECT(U = allocMatrix(REALSXP, nr, mindim));
  PROTECT(Vt = allocMatrix(REALSXP, mindim, nc));
  PROTECT(D = allocVector(REALSXP, mindim));

  /* duplicate the matrix so the original one is not overwritten. */
  if ((duplicated = NAMED(matrix)) > 0)
    PROTECT(matrix = duplicate(matrix));

  /* call Lapack to do the grunt work. */
  c_svd(REAL(matrix), REAL(U), REAL(D), REAL(Vt), &nr, &nc, &mindim);

  /* build the return value. */
  PROTECT(result = allocVector(VECSXP, 3));
  PROTECT(elnames = allocVector(STRSXP, 3));
  SET_STRING_ELT(elnames, 0, mkChar("d"));
  SET_STRING_ELT(elnames, 1, mkChar("u"));
  SET_STRING_ELT(elnames, 2, mkChar("vt"));
  setAttrib(result, R_NamesSymbol, elnames);
  SET_VECTOR_ELT(result, 0, D);
  SET_VECTOR_ELT(result, 1, U);
  SET_VECTOR_ELT(result, 2, Vt);

  if (duplicated > 0)
    UNPROTECT(6);
  else
    UNPROTECT(5);

  return result;

}/*R_SVD*/

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
static void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim) {

int err = 0, lwork = -1, *iwork = NULL;
char jobz = 'A';
double tmp = 0, *work = NULL;

  iwork = (int *)Calloc(8 * (*mindim), int);

  /* ask for the optimal size of the work array. */
  F77_CALL(dgesdd)(&jobz, nrows, ncols, A, nrows, D, U, nrows, 
                   V, mindim, &tmp, &lwork, iwork, &err);

  lwork = (int)tmp;
  work = (double *)Calloc(lwork, double);

  /* actual call */
  F77_NAME(dgesdd)(&jobz, nrows, ncols, A, nrows, D, U, nrows, 
                   V, mindim, work, &lwork, iwork, &err);

  Free(work);
  Free(iwork);

  if (err)
    error("an error (%d) occurred in the call to dgesdd().\n", err);

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

