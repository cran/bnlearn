#include "common.h"
#include <Rmath.h>
#include <R_ext/Lapack.h>

/* -------------------- some useful macros ------------------------------ */

#define NROWS(x) INTEGER(getAttrib(x, R_DimSymbol))[0]
#define NCOLS(x) INTEGER(getAttrib(x, R_DimSymbol))[1]

/* -------------------- function declarations --------------------------- */

static void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim);

/* -------------------- R level interfaces to LAPACK -------------------- */

/* Singular Value Decomposition (SVD). */
SEXP r_svd(SEXP matrix) {

  SEXP U, D, Vt, result, elnames;
  int nrows = 0, ncols = 0, mindim = 0;

  /* compute the dimensions of the SVD components. */
  nrows = NROWS(matrix);
  ncols = NCOLS(matrix);
  mindim = imin2(nrows, ncols);

  /* allocate the U, D and t(V) matrices. */
  PROTECT(U = allocMatrix(REALSXP, nrows, mindim));
  PROTECT(Vt = allocMatrix(REALSXP, mindim, ncols));
  PROTECT(D = allocVector(REALSXP, mindim));

  /* duplicate the matrix so the original one is not overwritten. */
  if (NAMED(matrix) > 0)
    PROTECT(matrix = duplicate(matrix));

  /* call Lapack to do the grunt work. */
  c_svd(REAL(matrix), REAL(U), REAL(D), REAL(Vt), &nrows, &ncols, &mindim);

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

  if (NAMED(matrix) > 0)
    UNPROTECT(6);
  else
    UNPROTECT(5);

  return result;

}/*R_SVD*/

/* -------------------- C level interfaces to LAPACK -------------------- */

/* C-level wrapper around the dgesvd() F77 routine. Note that the input
 * matrix A is overwritten by dgesvd(), so it's sensible to have a
 * backup copy in case it's needed later. */
static void c_svd(double *A, double *U, double *D, double *V, int *nrows,
    int *ncols, int *mindim) {

  int err, lwork = -1;
  char jobu, jobvt;
  double work1, *work;

  if (*nrows < *ncols) {

    jobu = 'A';
    jobvt = 'S';
		
  }/*THEN*/
  else {

    jobu = 'S';
    jobvt = 'A';

  }/*ELSE*/

  F77_CALL(dgesvd)(&jobu, &jobvt, nrows, ncols, A, nrows, D, U, nrows, V,
    mindim, &work1, &lwork, &err);

  lwork = (int)floor(work1);

  if (work1 - lwork > 0.5) lwork++;

  work = (double *)Calloc(lwork, double);

  /* actual call */
  F77_NAME(dgesvd)(&jobu, &jobvt, nrows, ncols, A, nrows, D, U, nrows, V,
    mindim, work, &lwork, &err);

  Free(work);

  if (err)
    error("an error (%d) occurred in the call to dgesvd().\n", err);

}/*C_SVD*/

