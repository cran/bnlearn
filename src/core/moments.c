#include "../include/rcore.h"
#include "allocations.h"
#include "moments.h"
#include "../minimal/common.h"

/* compute the sum of the squared errors. */
double c_sse(double *data, double mean, int nrow) {

int i = 0;
long double nvar = 0;

  for (i = 0; i < nrow; i++)
    nvar += (data[i] - mean) * (data[i] - mean);

  return (double)nvar;

}/*C_SSE*/

/* compute the mean. */
double c_mean(double *data, int nrow) {

int i = 0;
long double mean = 0;

  for (i = 0; i < nrow; i++)
    mean += data[i];
  mean /= nrow;

  return (double)mean;

}/*C_MEAN*/

/* compute a vector of means from a set of variables. */
void c_meanvec(double **data, double *mean, int nrow, int ncol, int first) {

int i = 0, j = 0;
long double sum = 0;

  for (i = first; i < ncol; i++) {

    for (j = 0, sum = 0; j < nrow; j++)
      sum += data[i][j];

    mean[i] = sum / nrow;

  }/*FOR*/

}/*C_MEANVEC*/

/* compute a vector of variances from a set of variables. */
void c_ssevec(double **data, double *sse, double *means, int nrow, int ncol,
    int first) {

int i = 0, j = 0;
long double sum = 0;

  for (i = first; i < ncol; i++) {

    for (sum = 0, j = 0 ; j < nrow; j++)
      sum += (data[i][j] - means[i]) * (data[i][j] - means[i]);

    sse[i] = sum;

  }/*FOR*/

}/*C_SSEVEC*/

/* compute a single standard deviation. */
void c_sd(double *xx, int nobs, int p, double mean, double *sd) {

  if (nobs == 0)
    *sd = NA_REAL;
  else if (nobs <= p)
    *sd = 0;
  else
    *sd = sqrt(c_sse(xx, mean, nobs) / (nobs - p));

}/*C_SD*/

/* compute one standard deviation per stratum. */
void c_cgsd(double *xx, int *z, int *nz, int nobs, int nstrata, int p,
  long double *means, double *sd) {

int i = 0, *nnz = NULL;
long double *ssr = NULL, *mm = NULL;

  /* allocate the accumulators. */
  ssr = Calloc1D(nstrata, sizeof(long double));

  /* allocate and compute the per-stratum sample sizes, if needed. */
  if (!nz) {

    nnz = Calloc1D(nstrata, sizeof(int));
    for (i = 0; i < nobs; i++)
      nnz[z[i] - 1]++;

  }/*THEN*/
  else {

    nnz = nz;

  }/*ELSE*/

  /* allocate and compute the per-stratum means, if needed. */
  if (!means) {

    mm = Calloc1D(nstrata, sizeof(long double));
    for (i = 0; i < nobs; i++)
      mm[z[i] - 1] += xx[i];
    for (i = 0; i < nstrata; i++)
      if (nnz[i] != 0)
        mm[i] /= nnz[i];

  }/*THEN*/
  else {

    mm = means;

  }/*ELSE*/

  /* compute the per-stratum sum of the squares residuals. */
  for (i = 0; i < nobs; i++)
    ssr[z[i] - 1] += (xx[i] - mm[z[i] - 1]) * (xx[i] - mm[z[i] - 1]);
  /* compute the standard deviations appropriately. */
  for (i = 0; i < nstrata; i++) {

    if (nnz[i] == 0)
      sd[i] = NA_REAL;
    else if (nnz[i] <= p)
      sd[i] = 0;
    else
      sd[i] = sqrt(ssr[i] / (nnz[i] - p));

  }/*FOR*/

  if (!nz)
    Free1D(nnz);
  if (!means)
    Free1D(mm);
  Free1D(ssr);

}/*C_CGSD*/

/* compute one or more standard deviations, possibly with strata. */
SEXP cgsd(SEXP x, SEXP strata, SEXP nparams) {

int nstrata = 0, nobs = length(x), *z = NULL;
double mean = 0, *xx = REAL(x);
SEXP sd;

  if (strata == R_NilValue) {

    PROTECT(sd = allocVector(REALSXP, 1));
    mean = c_mean(xx, nobs);
    c_sd(xx, nobs, INT(nparams), mean, REAL(sd));

  }/*THEN*/
  else {

    nstrata = NLEVELS(strata);
    z = INTEGER(strata);
    PROTECT(sd = allocVector(REALSXP, nstrata));

    c_cgsd(xx, z, NULL, nobs, nstrata, INT(nparams), NULL, REAL(sd));

  }/*ELSE*/

  UNPROTECT(1);

  return sd;

}/*CGSD*/

