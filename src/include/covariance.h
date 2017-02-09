
/* from covariance.c */
#define SD_GUARD(n, p, ret, expr) \
  if ((n) == 0) \
    (ret) = R_NaN; \
  else if ((n) <= (p)) \
    (ret) = 0; \
  else { \
    expr; \
  }/*SD_GUARD*/

double c_sse(double *data, double mean, int nrow);
double c_mean(double *data, int nrow);
void c_meanvec(double **data, double *mean, int nrow, int ncol, int first);
void c_ssevec(double **data, double *sse, double *means, int nrow, int ncol,
    int first);
void c_update_meanvec(double **data, double *mean, int update, int nrow);
void c_covmat(double **data, double *mean, int nrow, int ncol, double *mat,
    int first);
void c_update_covmat(double **data, double *mean, int update, int nrow,
    int ncol, double *mat);
void c_sd(double *xx, int nobs, int p, double mean, int compute, double *sd);
void c_cgsd(double *xx, int *z, int *nz, int nobs, int nstrata, int p,
    long double *means, double *sd);

/* from linear.correlation.c */
double c_fast_cor(double *xx, double *yy, int num, double xm, double ym,
    long double xsd, long double ysd);
double c_fast_pcor(double *cov, double *u, double *d, double *vt, int ncol,
    int strict);

