
/* from covariance.c */
double c_sse(double *data, double mean, int nrows);
double c_mean(double *data, int nrows);
void c_meanvec(double **data, double *mean, int nrows, int ncols, int first);
void c_ssevec(double **data, double *sse, double *means, int nrows, int ncols,
    int first);
void c_update_meanvec(double **data, double *mean, int update, int nrows);
void c_covmat(double **data, double *mean, int ncols, int nrows, double *mat,
    int first);
void c_update_covmat(double **data, double *mean, int update, int ncols,
    int nrows, double *mat);

/* from linear.correlation.c */
double c_fast_cor(double *xx, double *yy, int num, double xm, double ym,
    long double xsd, long double ysd);
double c_fast_pcor(double *cov, double *u, double *d, double *vt, int ncols,
    int strict);

