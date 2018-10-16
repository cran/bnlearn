#include "include/rcore.h"
#include "include/sampling.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/tests.h"
#include "include/blas.h"
#include "include/matrix.h"

/* unconditional Monte Carlo simulation for correlation-based tests. */
static double mc_cov(double *xx, double *yy, double xm, double ym, int n) {

int i = 0;
double sum = 0;

  /* compute the actual covariance. */
  for (i = 0; i < n; i++)
    sum += (xx[i] - xm) * (yy[i] - ym);

  return sum;

}/*MC_COV*/

void c_gauss_mcarlo(double *xx, double *yy, int num, int B, double *res,
    double alpha, test_e test, double *observed) {

int j = 0, k = 0;
double *yperm = NULL;
double enough = ceil(alpha * B) + 1, xm = 0, ym = 0, xsse = 0, ysse = 0;
int *perm = NULL, *work = NULL;

  /* cache the means of the two variables (invariant under permutation). */
  for (j = 0; j < num; j++) {

    xm += xx[j];
    ym += yy[j];

  }/*FOR*/

  xm /= num;
  ym /= num;

  /* compute the variances of the two variables (also invariant). */
  xsse = c_sse(xx, xm, num);
  ysse = c_sse(yy, ym, num);

  /* if at least one of the two variables is constant, they are independent. */
  if ((xsse == 0) || (ysse == 0)) {

    *observed = 0;
    *res = 1;

    return;

  }/*THEN*/

  /* allocate the arrays needed by RandomPermutation. */
  perm = Calloc1D(num, sizeof(int));
  work = Calloc1D(num, sizeof(int));
  /* allocate the array for the pemutations. */
  yperm = Calloc1D(num, sizeof(double));


  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random permutations of the second variable and check how many tests are
     greater (in absolute value) than the original one.*/
  *observed = mc_cov(xx, yy, xm, ym, num);

  for (j = 0; j < B; j++) {

    RandomPermutation(num, perm, work);

    for (k = 0; k < num; k++)
      yperm[k] = yy[perm[k] - 1];

    if (fabs(mc_cov(xx, yperm, xm, ym, num)) > fabs(*observed)) {

      SEQUENTIAL_COUNTER_CHECK(*res);

    }/*THEN*/

  }/*FOR*/

  /* compute the observed value for the statistic. */
  switch(test) {

    case MC_MI_G:
    case SMC_MI_G:
      *observed = c_fast_cor(xx, yy, num, xm, ym, xsse, ysse);
      *observed = 2 * num * cor_mi_trans(*observed);
      break;

    case MC_COR:
    case SMC_COR:
      *observed = c_fast_cor(xx, yy, num, xm, ym, xsse, ysse);
      break;

    case MC_ZF:
    case SMC_ZF:

      /* check whether the sample size is big enough for the transform. */
      if (num - 3 < 1) {

        warning("sample size too small to compute the Fisher's Z transform.");
        *observed = 0;

      }/*THEN*/
      else {

        *observed = cor_zf_trans(c_fast_cor(xx, yy, num, xm, ym, xsse, ysse),
                      (double)num - 2);

      }/*ELSE*/

      break;

    default:
      error("unknown permutation test statistic.");

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed p-value. */
  *res /= B;

  Free1D(perm);
  Free1D(work);
  Free1D(yperm);

}/*C_GAUSS_MCARLO*/

/* conditional Monte Carlo simulation for correlation-based tests. */
void c_gauss_cmcarlo(double **column, int ncol, int num, int v1, int v2, int B,
    double *observed, double *pvalue, double alpha, test_e test) {

int j = 0, k = 0, errcode = 0, *work = NULL, *perm = NULL;
int error_counter = 0;
double permuted = 0, *yperm = NULL, *yorig = NULL;
double enough = ceil(alpha * B) + 1;
double *mean = NULL;
covariance cov = { 0 }, backup = { 0 };

  /* cache the means of the variables (they are invariant under permutation). */
  mean = Calloc1D(ncol, sizeof(double));
  /* compute the mean values  */
  c_meanvec(column, mean, num, ncol, 0);

  /* allocate and initialize the covariance matrix. */
  cov = new_covariance(ncol, TRUE);
  backup = new_covariance(ncol, TRUE);
  c_covmat(column, mean, num, ncol, cov, 0);

  /* if at least one of the two variables is constant, they are independent. */
  if ((cov.mat[CMC(v1, v1, ncol)] == 0) || (cov.mat[CMC(v2, v2, ncol)] == 0)) {

    *observed = 0;
    *pvalue = 1;

    Free1D(mean);
    FreeCOV(backup);
    FreeCOV(cov);

    return;

  }/*THEN*/

  /* make a backup copy that will not be touched by permutations. */
  copy_covariance(&cov, &backup);

  /* substitute the original data with the fake column that will be permuted. */
  yperm = Calloc1D(num, sizeof(double));
  yorig = column[v2];
  memcpy(yperm, yorig, num * sizeof(double));
  column[v2] = yperm;

   /* allocate the arrays needed by RandomPermutation. */
  perm = Calloc1D(num, sizeof(int));
  work = Calloc1D(num, sizeof(int));

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random permutations (all variable but the second are fixed) and check how
     many tests are greater (in absolute value) than the original one.*/
  *observed = c_fast_pcor(cov, v1, v2, &errcode, TRUE);

  if (errcode)
    error("an error (%d) occurred in the call to dgesvd().\n", errcode);

  for (j = 0; j < B; j++) {

    /* reset the error flag of the SVD Fortran routine. */
    errcode = 0;

    RandomPermutation(num, perm, work);

    for (k = 0; k < num; k++)
      yperm[k] = yorig[perm[k] - 1];

    /* restore the covariance matrix from the good copy. */
    copy_covariance(&backup, &cov);
    /* update the relevant covariances. */
    c_update_covmat(column, mean, v2, num, ncol, cov.mat);

    permuted = c_fast_pcor(cov, v1, v2, &errcode, TRUE);

    if (errcode != 0)
      error_counter++;

    if (fabs(permuted) > fabs(*observed)) {

      SEQUENTIAL_COUNTER_CHECK(*pvalue);

    }/*THEN*/

  }/*FOR*/

  if (error_counter > 0)
    warning("unable to compute %d permutations due to errors in dgesvd().\n",
      error_counter);

  /* compute the observed value for the statistic. */
  switch(test) {

    case MC_MI_G:
    case SMC_MI_G:
      *observed = 2 * num * cor_mi_trans(*observed);
      break;

    case MC_COR:
    case SMC_COR:
      break;

    case MC_ZF:
    case SMC_ZF:
      /* check whether the sample size is big enough for the transform. */
      if (num - 1 - ncol < 1) {

        warning("sample size too small to compute the Fisher's Z transform.");
        *observed = 0;

      }/*THEN*/
      else {

        *observed = cor_zf_trans(*observed, (double)num - ncol);

      }/*ELSE*/

      break;

    default:
      error("unknown permutation test statistic.");

  }/*SWITCH*/

  PutRNGstate();

  /* restore the pointer to the original column. */
  column[v2] = yorig;

  /* save the observed p-value. */
  *pvalue /= B;

  Free1D(mean);
  Free1D(perm);
  Free1D(work);
  Free1D(yperm);
  FreeCOV(backup);
  FreeCOV(cov);

}/*C_GAUSS_CMCARLO*/

