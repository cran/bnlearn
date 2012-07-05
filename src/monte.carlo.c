
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "common.h"

#define MUTUAL_INFORMATION             1
#define PEARSON_X2                     2
#define GAUSSIAN_MUTUAL_INFORMATION    3
#define LINEAR_CORRELATION             4
#define FISHER_Z                       5

/* initialize the table of log-factorials. */
#define allocfact(n) \
  fact = alloc1dreal(n + 1); \
  fact[0] = 0.; \
  for(k = 1; k <= n; k++) \
    fact[k] = lgammafn((double) (k + 1));

#define sequential_counter_check(counter) \
  (counter)++; \
  if ((counter) >= enough) { \
    (counter) = *B; \
    break; \
  }

/* function declarations of the custom test functions. */
static double _mi(int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length);
static double _cmi(int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl);
static double _x2(int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length);
static double _cx2(int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl);
static double _cov(double *xx, double *yy, double *xm, double *ym, int *n);

/* unconditional Monte Carlo simulation for discrete tests. */
SEXP mcarlo(SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length, SEXP samples,
    SEXP test, SEXP alpha) {

double *fact = NULL, *res = NULL, observed = 0;
int *n = NULL, *ncolt = NULL, *nrowt = NULL, *workspace = NULL;
int *num = INTEGER(length), *nr = INTEGER(lx), *nc = INTEGER(ly);
int *xx = INTEGER(x), *yy = INTEGER(y), *B = INTEGER(samples);
int i = 0, k = 0, enough = ceil(NUM(alpha) * (*B)) + 1;
SEXP result;

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  res[0] = res[1] = 0;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(*num);

  /* allocate and initialize the workspace for rcont2. */
  workspace = alloc1dcont(*nc);

  /* initialize the contingency table. */
  n = alloc1dcont(*nr * (*nc));

  /* initialize the marginal frequencies. */
  nrowt = alloc1dcont(*nr);
  ncolt = alloc1dcont(*nc);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++)
    n[CMC(xx[k] - 1, yy[k] - 1, *nr)]++;

  /* compute the marginals. */
  for (i = 0; i < *nr; i++)
    for (k = 0; k < *nc; k++) {

      nrowt[i] += n[CMC(i, k, *nr)];
      ncolt[k] += n[CMC(i, k, *nr)];

    }/*FOR*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(INT(test)) {

    case MUTUAL_INFORMATION:
      observed = _mi(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < *B; k++) {

        rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (_mi(n, nrowt, ncolt, nr, nc, num) > observed) {

          sequential_counter_check(res[1]);

        }/*THEN*/

      }/*FOR*/

      observed = 2 * observed;

      break;

    case PEARSON_X2:
      observed = _x2(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < *B; k++) {

        rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (_x2(n, nrowt, ncolt, nr, nc, num) > observed) {

          sequential_counter_check(res[1]);

        }/*THEN*/

      }/*FOR*/

      break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed value of the statistic and the corresponding p-value. */
  NUM(result) =  observed;
  REAL(result)[1] =  REAL(result)[1] / (*B);

  UNPROTECT(1);

  return result;

}/*MCARLO*/

/* conditional Monte Carlo simulation for discrete tests. */
SEXP cmcarlo(SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz,
    SEXP length, SEXP samples, SEXP test, SEXP alpha) {

double *fact = NULL, *res = NULL, observed = 0;
int **n = NULL, **ncolt = NULL, **nrowt = NULL, *ncond = NULL, *workspace = NULL;
int *num = INTEGER(length), *B = INTEGER(samples);
int *nr = INTEGER(lx), *nc = INTEGER(ly), *nl = INTEGER(lz);
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
int i = 0, j = 0, k = 0, enough = ceil(NUM(alpha) * (*B)) + 1;
SEXP result;

  /* allocate and initialize the result */
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  res[0] = res[1] = 0;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(*num);

  /* allocate and initialize the workspace for rcont2. */
  workspace = alloc1dcont(*nc);

  /* initialize the contingency table. */
  n = alloc2dcont(*nl, (*nr) * (*nc));

  /* initialize the marginal frequencies. */
  nrowt = alloc2dcont(*nl, *nr);
  ncolt = alloc2dcont(*nl, *nc);
  ncond = alloc1dcont(*nl);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++)
    n[zz[k] - 1][CMC(xx[k] - 1, yy[k] - 1, *nr)]++;

  /* compute the marginals. */
  for (i = 0; i < *nr; i++)
    for (j = 0; j < *nc; j++)
      for (k = 0; k < *nl; k++) {

        nrowt[k][i] += n[k][CMC(i, j, *nr)];
        ncolt[k][j] += n[k][CMC(i, j, *nr)];
        ncond[k] += n[k][CMC(i, j, *nr)];

      }/*FOR*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(INT(test)) {

    case MUTUAL_INFORMATION:
      observed = _cmi(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < *B; j++) {

        for (k = 0; k < *nl; k++)
          rcont2(nr, nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (_cmi(n, nrowt, ncolt, ncond, nr, nc, nl) > observed) {

          sequential_counter_check(res[1]);

        }/*THEN*/


      }/*FOR*/

      observed = 2 * observed;

      break;

    case PEARSON_X2:
      observed = _cx2(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < *B; j++) {

        for (k = 0; k < *nl; k++)
          rcont2(nr, nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (_cx2(n, nrowt, ncolt, ncond, nr, nc, nl) > observed) {

          sequential_counter_check(res[1]);

        }/*THEN*/

      }/*FOR*/

      break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed value of the statistic and the corresponding p-value. */
  res[0] = observed;
  res[1] /= *B;

  UNPROTECT(1);

  return result;

}/*CMCARLO*/

/* unconditional Monte Carlo simulation for correlation-based tests. */
SEXP gauss_mcarlo(SEXP x, SEXP y, SEXP samples, SEXP test, SEXP alpha) {

int j = 0, k = 0, num = LENGTH(x), *B = INTEGER(samples);
double *xx = REAL(x), *yy = REAL(y), *yperm = NULL, *res = NULL;
double observed = 0, enough = ceil(NUM(alpha) * (*B)) + 1, xm = 0, ym = 0;
int *perm = NULL, *work = NULL;
SEXP result;

  /* allocate the arrays needed by RandomPermutation. */
  perm = alloc1dcont(num);
  work = alloc1dcont(num);

  /* allocate the array for the pemutations. */
  yperm = alloc1dreal(num);

  /* cache the means of the two variables (they are invariant under permutation). */
  for (j = 0; j < num; j++) {

    xm += xx[j];
    ym += yy[j];

  }/*FOR*/

  xm /= num;
  ym /= num;

  /* allocate the result. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random permutations (all variable but the second are fixed) and check how
     many tests are greater (in absolute value) than the original one.*/
  switch(INT(test)) {

    case GAUSSIAN_MUTUAL_INFORMATION:
    case LINEAR_CORRELATION:
    case FISHER_Z:
      observed = _cov(xx, yy, &xm, &ym, &num);

      for (j = 0; j < *B; j++) {

        RandomPermutation(num, perm, work);

        for (k = 0; k < num; k++)
          yperm[k] = yy[perm[k] - 1];

        if (fabs(_cov(xx, yperm, &xm, &ym, &num)) > fabs(observed)) {

          sequential_counter_check(*res);

        }/*THEN*/

      }/*FOR*/

    break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed p-value. */
  *res /= *B;

  UNPROTECT(1);

  return result;

}/*GAUSS_MCARLO*/

/* conditional Monte Carlo simulation for correlation-based tests. */
SEXP gauss_cmcarlo(SEXP data, SEXP length, SEXP samples, SEXP test, SEXP alpha) {

int j = 0, k = 0, ncols = LENGTH(data), errcode = 0, *work = NULL, *perm = NULL;
int error_counter = 0, *B = INTEGER(samples), *num = INTEGER(length);
double observed = 0, permuted = 0, *yperm = NULL, *yorig = NULL, *res = NULL;
double enough = ceil(NUM(alpha) * (*B)) + 1;
double **column = NULL, *mean = NULL, *covariance = NULL, *covariance_backup = NULL;
double *u = NULL, *d = NULL, *vt = NULL;
SEXP result;

  /* allocate the matrices needed for the SVD decomposition. */
  u = alloc1dreal(ncols * ncols);
  d = alloc1dreal(ncols);
  vt = alloc1dreal(ncols * ncols);

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* allocate and initialize an array of pointers for the variables. */
  column = (double **) alloc1dpointer(ncols);
  for (j = 0; j < ncols; j++)
    column[j] = REAL(VECTOR_ELT(data, j));

  /* cache the means of the variables (they are invariant under permutation). */
  mean = alloc1dreal(ncols);

  /* compute the mean values  */
  for (j = 0; j < ncols; j++) {

    for (k = 0 ; k < *num; k++)
      mean[j] += column[j][k];

    mean[j] /= (*num);

  }/*FOR*/

  /* allocate and initialize the covariance matrix. */
  covariance = alloc1dreal(ncols * ncols);
  covariance_backup = alloc1dreal(ncols * ncols);
  c_covmat(column, mean, &ncols, num, covariance);
  memcpy(covariance_backup, covariance, ncols * ncols * sizeof(double));

  /* substitute the original data with the fake column that will be permuted. */
  yperm = alloc1dreal(*num);
  yorig = column[1];
  memcpy(yperm, yorig, *num * sizeof(double));
  column[1] = yperm;

   /* allocate the arrays needed by RandomPermutation. */
  perm = alloc1dcont(*num);
  work = alloc1dcont(*num);

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random permutations (all variable but the second are fixed) and check how
     many tests are greater (in absolute value) than the original one.*/
  switch(INT(test)) {

    case GAUSSIAN_MUTUAL_INFORMATION:
    case LINEAR_CORRELATION:
    case FISHER_Z:
      observed = c_fast_pcor(covariance, &ncols, u, d, vt, &errcode);

      if (errcode)
        error("an error (%d) occurred in the call to dgesvd().\n", errcode);

      for (j = 0; j < (*B); j++) {

        /* reset the error flag of the SVD Fortran routine. */
        errcode = 0;

        RandomPermutation(*num, perm, work);

        for (k = 0; k < *num; k++)
          yperm[k] = yorig[perm[k] - 1];

        /* restore the covariance matrix from the good copy. */
        memcpy(covariance, covariance_backup, ncols * ncols * sizeof(double));
        /* update the relevant covariances. */
        c_update_covmat(column, mean, 1, &ncols, num, covariance);

        permuted = c_fast_pcor(covariance, &ncols, u, d, vt, &errcode);

        if (errcode != 0)
          error_counter++;

        if (fabs(permuted) > fabs(observed)) {

          sequential_counter_check(*res);

        }/*THEN*/

      }/*FOR*/

    if (error_counter > 0)
      warning("unable to compute %d permutations due to errors in dgesvd().\n",
        error_counter);

    break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed p-value. */
  *res /= *B;

  UNPROTECT(1);

  return result;

}/*GAUSS_CMCARLO*/

/* compute the mutual information from the joint and marginal frequencies. */
static double _mi(int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < *nrows; i++)
    for (j = 0; j < *ncols; j++)
      res += MI_PART(n[CMC(i, j, *nrows)], nrowt[i], ncolt[j], *length); 

  return res;

}/*_MI*/

/* compute the conditional mutual information from the joint and marginal frequencies. */
static double _cmi(int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < *nl; k++)
    for (j = 0; j < *nc; j++)
      for (i = 0; i < *nr; i++) 
        res += MI_PART(n[k][CMC(i, j, *nr)], nrowt[k][i], ncolt[k][j], ncond[k]);

  return res;

}/*_CMI*/

/* compute Pearson's X^2 coefficient from the joint and marginal frequencies. */
static double _x2(int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < *nrows; i++)
    for (j = 0; j < *ncols; j++) {

      if (n[CMC(i, j, *nrows)] != 0)
        res += (n[CMC(i, j, *nrows)] - nrowt[i] * (double)ncolt[j] / (*length)) *
               (n[CMC(i, j, *nrows)] - nrowt[i] * (double)ncolt[j] / (*length)) /
               (nrowt[i] * (double)ncolt[j] / (*length));

    }/*FOR*/

  return res;

}/*_X2*/

/* compute the Pearson's conditional X^2 coefficient from the joint and marginal frequencies. */
static double _cx2(int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < *nl; k++)
    for (j = 0; j < *nc; j++)
      for (i = 0; i < *nr; i++) {

       if (n[k][CMC(i, j, *nr)] != 0) {

          res += (n[k][CMC(i, j, *nr)] - nrowt[k][i] * (double)ncolt[k][j] / ncond[k]) *
                 (n[k][CMC(i, j, *nr)] - nrowt[k][i] * (double)ncolt[k][j] / ncond[k]) /
                 (nrowt[k][i] * (double)ncolt[k][j] / ncond[k]);

        }/*THEN*/

      }/*FOR*/

  return res;

}/*_CX2*/

/* compute a (barebone version of) the linear correlation coefficient. */
static double _cov(double *xx, double *yy, double *xm, double *ym, int *n) {

int i = 0;
double sum = 0;

  /* compute the actual covariance. */
  for (i = 0; i < *n; i++)
    sum += (xx[i] - *xm) * (yy[i] - *ym);

  return sum;

}/*_COV*/

