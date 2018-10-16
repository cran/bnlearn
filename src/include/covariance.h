#include "data.structures.h"

/* from covariance.c */
covariance new_covariance(int dim, int decomp);
void copy_covariance(covariance *src, covariance *copy);
void covariance_drop_column(covariance *cov, int to_drop);
void FreeCOV(covariance cov);

double c_sse(double *data, double mean, int nrow);
double c_mean(double *data, int nrow);
void c_meanvec(double **data, double *mean, int nrow, int ncol, int first);
void c_ssevec(double **data, double *sse, double *means, int nrow, int ncol,
    int first);
void c_covmat(double **data, double *mean, int nrow, int ncol, covariance cov,
    int first);
void c_update_covmat(double **data, double *mean, int update, int nrow,
    int ncol, double *mat);
void c_covmat_with_missing(double **data, int nrow, int ncol,
    short int *missing_yz, short int *missing_all, double *mean, double *mat,
    int *ncomplete);
void c_sd(double *xx, int nobs, int p, double mean, int compute, double *sd);
void c_cgsd(double *xx, int *z, int *nz, int nobs, int nstrata, int p,
    long double *means, double *sd);

/* from linear.correlation.c */
double c_fast_cor(double *xx, double *yy, int num, double xm, double ym,
    long double xsd, long double ysd);
double c_cor_with_missing(double *x, double *y, int nobs, double *xm,
    double *ym, double *xsd, double *ysd, int *ncomplete);
double c_fast_pcor(covariance cov, int v1, int v2, int *err, int decomp);
