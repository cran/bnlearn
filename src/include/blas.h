
/* from linear.algebra.c */
SEXP r_svd(SEXP matrix, SEXP strict);
SEXP r_det(SEXP matrix, int scale);
double c_det(double *matrix, int *rows);
void c_udvt(double **u, double **d, double **vt, int ncol);
void c_svd(double *A, double *U, double *D, double *V, int *nrow, int *ncol,
    int *mindim, int strict, int *errcode);
void c_ginv(double *covariance, int ncol, double *mpinv);
void c_finv(double *cov, int *ncol, double *mpinv);
double c_quadratic(double *x, int *ncol, double *sigma, double *y,
    double *workspace);
void c_rotate(double *S1, double *S2, double *x, double *a, double *mu,
    int *ncol, double *workspace);
void c_qr(double *qr, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd);

/* from least.squares.c */
void c_ols(double **x, double *y, int nrow, int ncol, double *fitted,
    double *resid, double *beta, double *sd, int missing);
void c_cls(double **x, double *y, int *z, int nrow, int ncol, int ncond,
    double *fitted, double *resid, double *beta, double *sd, int missing);

/* from common.c */
void c_qr_matrix(double *qr, double **x, int nrow, int ncol, int *complete,
    int ncomplete);

