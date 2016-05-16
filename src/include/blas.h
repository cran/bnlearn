
/* from linear.algebra.c */
SEXP r_svd(SEXP matrix, SEXP strict);
SEXP r_det(SEXP matrix, int scale);
double c_det(double *matrix, int *rows);
void c_udvt(double **u, double **d, double **vt, int ncols);
void c_svd(double *A, double *U, double *D, double *V, int *nrows, int *ncols,
    int *mindim, int strict, int *errcode);
void c_ginv(double *covariance, int ncols, double *mpinv);
void c_finv(double *cov, int *ncols, double *mpinv);
double c_quadratic(double *x, int *ncols, double *sigma, double *y,
    double *workspace);
void c_rotate(double *S1, double *S2, double *x, double *a, double *mu,
    int *ncols, double *workspace);
void c_qr_ols (double *qr, double *y, int nrow, int ncol, double *fitted,
    long double *sd);

/* from common.c */
SEXP qr_matrix(SEXP dataframe, SEXP name);

