
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
  fact = alloc1dreal(n); \
  fact[0] = 0.; \
  for(k = 1; k <= n; k++) \
    fact[k] = lgammafn((double) (k + 1)); \

/* function declarations of the custom test functions. */
static double _mi (int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length);
static double _cmi (int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl);
static double _x2 (int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length);
static double _cx2 (int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl);
static double _cov(double *xx, double *yy, int *n);

/* unconditional Monte Carlo simulation for discrete tests. */
SEXP mcarlo (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length, SEXP samples,
    SEXP test) {

  double *fact, observed = 0;
  int *n, *ncolt, *nrowt, *workspace;
  int *num = INTEGER(length), *nr = INTEGER(lx), *nc = INTEGER(ly);
  int *xx = INTEGER(x), *yy = INTEGER(y);
  int *B = INTEGER(samples);
  int k = 0;

  SEXP result;

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(REALSXP, 2));
  REAL(result)[1] = 0.;

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
  for (k = 0; k < *num; k++) {

    n[CMC(xx[k] - 1, yy[k] - 1, *nr)]++;
    nrowt[xx[k] - 1]++;
    ncolt[yy[k] - 1]++;

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

        if (_mi(n, nrowt, ncolt, nr, nc, num) > observed)
          REAL(result)[1] += 1;

      }/*FOR*/

      observed = 2 * observed;

      break;

    case PEARSON_X2:
      observed = _x2(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < *B; k++) {

        rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (_x2(n, nrowt, ncolt, nr, nc, num) > observed)
          REAL(result)[1] += 1;

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
SEXP cmcarlo (SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz,
    SEXP length, SEXP samples, SEXP test) {

  double *fact, *res, observed = 0;
  int **n, **ncolt, **nrowt, *ncond, *workspace;
  int *num = INTEGER(length), *B = INTEGER(samples);;
  int *nr = INTEGER(lx), *nc = INTEGER(ly), *nl = INTEGER(lz);
  int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
  int j = 0, k = 0;

  SEXP result;

  /* allocate and initialize the result */
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  res[1] = 0;

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
  for (k = 0; k < *num; k++) {

    n[zz[k] - 1][CMC(xx[k] - 1, yy[k] - 1, *nr)]++;
    nrowt[zz[k] - 1][xx[k] - 1]++;
    ncolt[zz[k] - 1][yy[k] - 1]++;
    ncond[zz[k] - 1]++;

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

        if (_cmi(n, nrowt, ncolt, ncond, nr, nc, nl) > observed)
          res[1] += 1;

      }/*FOR*/

      observed = 2 * observed;

      break;

    case PEARSON_X2:
      observed = _cx2(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < *B; j++) {

        for (k = 0; k < *nl; k++)
          rcont2(nr, nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (_cx2(n, nrowt, ncolt, ncond, nr, nc, nl) > observed)
          res[1] += 1;

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
SEXP gauss_mcarlo (SEXP x, SEXP y, SEXP samples, SEXP test) {

  int j = 0, k = 0;
  double *xx = REAL(x), *yy = REAL(y), *yperm, *res;
  int num = LENGTH(x), *B = INTEGER(samples);
  double observed = 0;
  int *perm, *work;
  SEXP result;

  /* allocate the arrays needed by RandomPermutation. */
  perm = alloc1dcont(num);
  work = alloc1dcont(num);

  /* allocate the array for the pemutations. */
  yperm = alloc1dreal(num);

  /* allocate the result. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater (in absolute value) than the original one.*/
  switch(INT(test)) {

    case GAUSSIAN_MUTUAL_INFORMATION:
    case LINEAR_CORRELATION:
    case FISHER_Z:
      observed = _cov(xx, yy, &num);

      for (j = 0; j < *B; j++) {

        RandomPermutation(num, perm, work);

        for (k = 0; k < num; k++)
          yperm[k] = yy[perm[k]];

        if (fabs(_cov(xx, yperm, &num)) > fabs(observed))
          *res += 1;

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
SEXP gauss_cmcarlo (SEXP data, SEXP length, SEXP samples, SEXP test) {

  SEXP data2, yy, result;
  int j = 0, k = 0;
  int *work, *perm;
  int *B = INTEGER(samples), *num = INTEGER(length);
  double observed = 0;
  double *yperm, *yorig, *res;

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* create a fake dataset which references all the columns of the original
   * data except the second one (which created from scratch and changed at
   * each iteration). */
  PROTECT(data2 = allocVector(VECSXP, LENGTH(data)));
  PROTECT(yy = allocVector(REALSXP, *num));
  yperm = REAL(yy);
  yorig = REAL(VECTOR_ELT(data, 1));

  for (j = 0; j < LENGTH(data); j++)
    if (j != 1)
      SET_VECTOR_ELT(data2, j, VECTOR_ELT(data, j));
    else
      SET_VECTOR_ELT(data2, j, yy);

   /* allocate the arrays needed by RandomPermutation. */
  perm = alloc1dcont(*num);
  work = alloc1dcont(*num);

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater (in absolute value) than the original one.*/
  switch(INT(test)) {

    case GAUSSIAN_MUTUAL_INFORMATION:
    case LINEAR_CORRELATION:
    case FISHER_Z:
      observed = NUM(fast_pcor(data, length));

      for (j = 0; j < (*B); j++) {

        RandomPermutation(*num, perm, work);
        
        for (k = 0; k < *num; k++)
          yperm[k] = yorig[perm[k]];

        if (fabs(NUM(fast_pcor(data2, length))) > fabs(observed))
          *res += 1;
 
      }/*FOR*/

    break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the observed p-value. */
  *res /= *B;

  UNPROTECT(3);

  return result;

}/*GAUSS_CMCARLO*/

/* compute the mutual information from the joint and marginal frequencies. */
static double _mi (int *n, int *nrowt, int *ncolt, int *nrows,
    int *ncols, int *length) {

  int i = 0, j = 0;
  double res = 0;

  for (i = 0; i < *nrows; i++)
    for (j = 0; j < *ncols; j++) {

      if (n[CMC(i, j, *nrows)] != 0)
        res += ((double)n[CMC(i, j, *nrows)]) *
                log((double)n[CMC(i, j, *nrows)]*(*length)/(double)(nrowt[i]*ncolt[j]));

    }/*FOR*/

  return res;

}/*_MI*/

/* compute the conditional mutual information from the joint and marginal frequencies. */
static double _cmi (int **n, int **nrowt, int **ncolt, int *ncond,
    int *nr, int *nc, int *nl) {

  int i = 0, j = 0, k = 0;
  double res = 0;

  for (k = 0; k < *nl; k++)
    for (j = 0; j < *nc; j++)
      for (i = 0; i < *nr; i++) {

       if (n[k][CMC(i, j, *nr)] != 0) {

          res += (double)n[k][CMC(i, j, *nr)] *
            log( (double)(n[k][CMC(i, j, *nr)]*ncond[k]) / (double)(nrowt[k][i] * ncolt[k][j]) );

        }/*THEN*/

      }/*FOR*/

  return res;

}/*_CMI*/

/* compute Pearson's X^2 coefficient from the joint and marginal frequencies. */
static double _x2 (int *n, int *nrowt, int *ncolt, int *nrows,
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
static double _cx2 (int **n, int **nrowt, int **ncolt, int *ncond,
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
static double _cov(double *xx, double *yy, int *n) {

  int i = 0;
  double sum = 0, xm = 0, ym = 0;

  /* compute the mean values  */
  for (i = 0; i < *n; i++) {

    xm += xx[i];
    ym += yy[i];

  }/*FOR*/

  xm /= *n;
  ym /= *n;

  /* compute the actual covariance. */
  for (i = 0; i < *n; i++)
    sum += (xx[i] - xm) * (yy[i] - ym);

  return sum;

}/*_COV*/

