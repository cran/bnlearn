#ifndef CORRELATION_HEADER
#define CORRELATION_HEADER

#include "covariance.matrix.h"

double c_fast_cor(double *xx, double *yy, int num, double xm, double ym,
    long double xsd, long double ysd);
double c_cor_with_missing(double *x, double *y, int nobs, double *xm,
    double *ym, double *xsd, double *ysd, int *ncomplete);
double c_fast_pcor(covariance cov, int v1, int v2, int *err, bool decomp);

#endif
