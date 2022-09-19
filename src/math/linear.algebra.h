#ifndef LINEAR_ALGEBRA_HEADER
#define LINEAR_ALGEBRA_HEADER

#include "../core/covariance.matrix.h"

/* column-major coordinates conversion macro. */
#define CMC(i, j, nrow) ((i) + (j) * (nrow))

/* matrix with elements stored in column-major order. */
typedef struct{

  int nrows;      /* first dimension. */
  int ncols;      /* second dimension */
  int *el;        /* pointer to the column-major matrix. */

} cmcmap;

#define CMEL(matrix, i, j) matrix.el[CMC(i, j, matrix.nrows)]

double c_logdet(double *matrix, int rows);
void c_udvt(double **u, double **d, double **vt, int ncol);
void c_svd(double *A, double *U, double *D, double *V, int *nrow, int *ncol,
    int *mindim, bool strict, int *errcode);
void c_ginv(covariance cov, covariance mpinv);
void c_qr(double *qr, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd);

void c_ols(double **x, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd, bool missing);
void c_cls(double **x, double *y, int *z, int nrow, int ncol, int ncond,
    double *fitted, double *resid, double *beta, double *sd, bool missing);

void c_qr_matrix(double *qr, double **x, int nrow, int ncol, int *complete,
    int ncomplete);

#endif
