#include "include/rcore.h"
#include "include/allocations.h"
#include "include/tests.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/matrix.h"

/* shrinked mutual information, to be used in C code. */
double c_shmi(int *xx, int llx, int *yy, int lly, int num) {

int i = 0, j = 0, k = 0;
double **n = NULL, *ni = NULL, *nj = NULL;
double lambda = 0, target = 1/(double)(llx * lly);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dreal(llx, lly);
  ni = alloc1dreal(llx);
  nj = alloc1dreal(lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;

  /* estimate the optimal lambda for the data. */
  mi_lambda((double *)n, &lambda, target, num, llx, lly, 0);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
        n[i][j] = lambda * target + (1 - lambda) * n[i][j] / num;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0)
        res += n[i][j] * log(n[i][j] / (ni[i] * nj[j]));

  return res;

}/*C_SHMI*/

/* shrinked conditional mutual information, to be used in C code. */
double c_shcmi(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df) {

int i = 0, j = 0, k = 0;
double ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
double lambda = 0, target = 1/(double)(llx * lly * llz);
double res = 0;

  /* compute the degrees of freedom. */
  *df = (double)(llx - 1) * (double)(lly - 1) * (double)(llz);

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dreal(llx, lly, llz);
  ni = alloc2dreal(llx, llz);
  nj = alloc2dreal(lly, llz);
  nk = alloc1dreal(llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;

  /* estimate the optimal lambda for the data. */
  mi_lambda((double *)n, &lambda, target, num, llx, lly, llz);

  /* switch to the probability scale and shrink the estimates. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++)
        n[i][j][k] = lambda * target + (1 - lambda) * n[i][j][k] / num;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  for (k = 0; k < llz; k++) {

    /* check each level of the conditioning variable to avoid (again)
     * "divide by zero" errors. */
    if (nk[k] == 0)
      continue;

    for (j = 0; j < lly; j++) {

      for (i = 0; i < llx; i++) {

        if (n[i][j][k] > 0)
          res += n[i][j][k] * log( (n[i][j][k] * nk[k]) / (ni[i][k] * nj[j][k]) );

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  return res;

}/*C_SHCMI*/

/* compute the shrinkage intensity lambda for the mutual information. */
void mi_lambda(double *n, double *lambda, double target, int num, int llx,
    int lly, int llz) {

double lden = 0, lnum = 0, temp = 0;

  /* compute the numerator and the denominator of the shrinkage intensity;
   * if the third dimension is a NULL pointer it's a 2-dimensional table. */
  if (llz == 0) {

    for (int i = 0; i < llx; i++)
      for (int j = 0; j < lly; j++) {

        temp = ((double **)n)[i][j] / (double)(num);
        lnum += temp * temp;
        temp = target - ((double **)n)[i][j] / (double)(num);
        lden += temp * temp;

      }/*FOR*/

  }/*THEN*/
  else {

    for (int i = 0; i < llx; i++)
      for (int j = 0; j < lly; j++)
        for (int k = 0; k < llz; k++) {

          temp = ((double ***)n)[i][j][k] / (double)(num);
          lnum += temp * temp;
          temp = target - ((double ***)n)[i][j][k] / (double)(num);
          lden += temp * temp;

      }/*FOR*/

  }/*ELSE*/

   /* compute the shrinkage intensity (avoiding "divide by zero" errors). */
  if (lden == 0)
    *lambda = 1;
  else
    *lambda = (1 - lnum) / ((double)(num - 1) * lden);

  /* bound the shrinkage intensity in the [0,1] interval. */
  if (*lambda > 1)
    *lambda = 1;
  if (*lambda < 0)
    *lambda = 0;

}/*MI_LAMBDA*/

/* shrinked linear correlation, to be used in C code. */
double c_fast_shcor(double *xx, double *yy, int *n) {

int i = 0;
double sum = 0;
double xm = 0, ym = 0, xsd = 0, ysd = 0, lambda = 0;
double tol = MACHINE_TOL;

  /* compute the mean values.  */
  for (i = 0 ; i < *n; i++) {

    xm += xx[i];
    ym += yy[i];

  }/*FOR*/

  xm /= (*n);
  ym /= (*n);

  /* compute the actual covariance. */
  for (i = 0; i < *n; i++) {

    sum += (xx[i] - xm) * (yy[i] - ym);
    xsd += (xx[i] - xm) * (xx[i] - xm);
    ysd += (yy[i] - ym) * (yy[i] - ym);

  }/*FOR*/

  /* note that the shrinkage intesity for the correlation coefficient is
   * identical to that for the covariance; so we don't need to standardize
   * the data. */
  for (i = 0; i < *n; i++) {

    lambda += ((xx[i] - xm) * (yy[i] - ym) - sum / (*n)) *
              ((xx[i] - xm) * (yy[i] - ym) - sum / (*n));

  }/*FOR*/

  /* compute lambda, the shrinkage intensity. */
  lambda *= (*n / sum) * (*n / sum);
  lambda *= (double)(*n) / (double)(*n - 1) / (double)(*n -1) / (double)(*n -1);

  /* truncate the shrinkage intensity in the [0,1] interval; this is not an
   * error, but a measure to increase the quality of the shrinked estimate. */
  if (lambda > 1) {

    lambda = 1;

  }/*THEN*/
  else if (lambda < 0) {

    lambda = 0;

  }/*THEN*/

  /* safety check against "divide by zero" errors. */
  if (fabs(xsd) < tol || fabs(ysd) < tol)
    sum = 0;
  else
    sum *= (1 - lambda) / (sqrt(xsd) * sqrt(ysd));

  /* double check that the coefficient is in the [-1, 1] range. */
  if (sum > 1) {

    warning("fixed correlation coefficient greater than 1, probably due to floating point errors.");

    sum = 1;

  }/*THEN*/
  else if (sum < -1) {

    warning("fixed correlation coefficient lesser than -1, probably due to floating point errors.");

    sum = -1;

  }/*ELSE*/

  return sum;

}/*C_FAST_SHCOR*/

/* Shrinked Covariance Matrix. */
void c_cov_lambda(double **column, double *mean, int ncols, int n, double *var) {

int i = 0, j = 0, k = 0, cur = 0;
double lambda = 0, sumcors = 0, sumvars = 0;

  for (i = 0; i < ncols; i++) {

    for (j = i; j < ncols; j++) {

      cur = CMC(i, j, ncols);

      /* compute the actual variance/covariance. */
      for (k = 0; k < n; k++)
        var[cur] += (column[i][k] - mean[i]) * (column[j][k] - mean[j]);

      if (i != j) {

        /* do the first round of computations for the shrinkage intensity. */
        for (k = 0; k < n; k++) {

          sumvars +=
            ((column[i][k] - mean[i]) * (column[j][k] - mean[j]) - var[cur] / n) *
            ((column[i][k] - mean[i]) * (column[j][k] - mean[j]) - var[cur] / n);

        }/*FOR*/

        sumcors += (var[cur] / (n - 1)) * (var[cur] / (n - 1));

      }/*THEN*/

      /* use the unbiased estimator for variances/covariances. */
      var[cur] /= n - 1;

      /* fill in the symmetric element of the matrix. */
      var[CMC(j, i, ncols)] = var[cur];

    }/*FOR*/

  }/*FOR*/

  /* wrap up the computation of the shrinkage intensity. */
  lambda = sumvars * n / (n - 1) / (n -1) / (n -1) / sumcors;

  /* truncate the shrinkage intensity in the [0,1] interval; this is not an
   * error, but a measure to increase the quality of the shrinked estimate. */
  if (lambda > 1) {

    lambda = 1;

  }/*THEN*/
  else if (lambda < 0) {

    lambda = 0;

  }/*THEN*/

  /* shrink the covariance matrix (except the diagonal, which stays the same). */
  for (i = 0; i < ncols; i++)
    for (j = 0; j < ncols; j++)
      if (i != j)
        var[CMC(i, j, ncols)] *= 1 - lambda;

}/*C_COV_LAMBDA*/

