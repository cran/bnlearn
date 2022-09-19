#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../tests.h"
#include "../../include/globals.h"
#include "../../core/covariance.matrix.h"
#include "../../math/linear.algebra.h"

#define TRUNCATE_LAMBDA(lambda) \
  if (lambda > 1) \
    lambda = 1; \
  if (lambda < 0) \
    lambda = 0;

/* compute the shrinkage intensity lambda for marginal correlation. */
double cor_lambda(double *xx, double *yy, int nobs, int ncomplete, double xm,
   double ym, double xsd, double ysd, double cor) {

int i = 0;
long double sum = 0, lambda = 0;

  sum = cor * sqrt(xsd * ysd) / (ncomplete - 1);

  /* note that the shrinkage intesity for the correlation coefficient is
   * identical to that for the covariance; so we don't need to standardize
   * the data. */
  for (i = 0; i < nobs; i++) {

    if (ISNAN(xx[i]) || ISNAN(yy[i]))
      continue;

    lambda += ((xx[i] - xm) * (yy[i] - ym) - sum) *
              ((xx[i] - xm) * (yy[i] - ym) - sum);

  }/*FOR*/

  if (lambda > MACHINE_TOL) {

    /* compute lambda, the shrinkage intensity, on a log-scale for numerical
     * stability (if lambda is equal to zero, just keep it as it is). */
    lambda = exp(log(lambda) - log(sum * sum) + log((double)ncomplete)
               - 3 * log((double)(ncomplete - 1)));

    /* truncate the shrinkage intensity in the [0,1] interval; this is not an
     * error, but a measure to increase the quality of the shrinked estimate. */
    TRUNCATE_LAMBDA(lambda);

  }/*THEN*/
  else {

    lambda = 0;

  }/*ELSE*/

  return (double)lambda;

}/*COR_LAMBDA*/

/* compute the shrinkage intensity lambda for a covariance matrix. */
double covmat_lambda(double **column, double *mean, covariance cov, int n,
    bool *missing, int nc) {

int i = 0, j = 0, k = 0, cur = 0;
long double lambda = 0, sum_covs = 0, sum_cov_vars = 0, temp = 0;

  for (i = 0; i < cov.dim; i++) {

    for (j = i; j < cov.dim; j++) {

      cur = CMC(i, j, cov.dim);

      /* only shrink off-diagonal elements. */
      if (i == j)
        continue;

      /* do the first round of computations for the shrinkage intensity. */
      for (k = 0; k < n; k++) {

        if (missing)
          if (missing[k])
            continue;

        temp = (column[i][k] - mean[i]) * (column[j][k] - mean[j]) -
                  (cov.mat[cur] * (double)(nc - 1) / (double)nc);
        sum_cov_vars += temp * temp;

      }/*FOR*/

      sum_covs += cov.mat[cur] * cov.mat[cur];

    }/*FOR*/

  }/*FOR*/

  if (sum_cov_vars > MACHINE_TOL) {

    /* compute lambda, the shrinkage intensity, on a log-scale for numerical
     * stability (if lambda is equal to zero, just keep it as it is). */
    lambda = exp(log(sum_cov_vars) + log((double)nc) - 3 * log((double)(nc - 1))
               - log(sum_covs));

    /* truncate the shrinkage intensity in the [0,1] interval; this is not an
     * error, but a measure to increase the quality of the shrinked estimate. */
    TRUNCATE_LAMBDA(lambda);

  }/*THEN*/
  else {

    lambda = 0;

  }/*ELSE*/

  return (double)lambda;

}/*COVMAT_LAMBDA*/

/* shrink the covariance matrix (except the diagonal, which stays the same). */
void covmat_shrink(covariance cov, double lambda) {

  for (int i = 0; i < cov.dim; i++)
    for (int j = 0; j < cov.dim; j++)
      if (i != j)
        cov.mat[CMC(i, j, cov.dim)] *= 1 - lambda;

}/*COVMAT_SHRINK*/

