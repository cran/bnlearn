#ifndef MOMENTS_HEADER
#define MOMENTS_HEADER

double c_sse(double *data, double mean, int nrow);
double c_mean(double *data, int nrow);

void c_meanvec(double **data, double *mean, int nrow, int ncol, int first);
void c_ssevec(double **data, double *sse, double *means, int nrow, int ncol,
    int first);

void c_sd(double *xx, int nobs, int p, double mean, double *sd);
void c_cgsd(double *xx, int *z, int *nz, int nobs, int nstrata, int p,
    long double *means, double *sd);

#endif
