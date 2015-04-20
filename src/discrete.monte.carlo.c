#include "include/rcore.h"
#include "include/allocations.h"
#include "include/sampling.h"
#include "include/tests.h"
#include "include/matrix.h"

/* initialize the table of log-factorials. */
#define allocfact(n) \
  fact = alloc1dreal(n + 1); \
  fact[0] = 0.; \
  for(k = 1; k <= n; k++) \
    fact[k] = lgammafn((double) (k + 1));

static double mc_mi(int *n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length);
static double mc_cmi(int **n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl);
static double mc_x2(int *n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length);
static double mc_cx2(int **n, int **nrowt, int **ncolt, int *ncond,
    int nr, int nc, int nl);
static double mc_jt(int *n, int *nrowt, int nrows, int ncols, int length);
static double mc_cjt(int **n, int **nrowt, int *ncond, int nr, int nc, int nl);
static double mc_cvjt(int **nrowt, int **ncolt, int *ncond, int nr, int nc,
    int nl);

/* unconditional Monte Carlo and semiparametric discrete tests. */
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, int test, double *df) {

double *fact = NULL;
int *n = NULL, *ncolt = NULL, *nrowt = NULL, *workspace = NULL;
int i = 0, k = 0, enough = ceil(alpha * B) + 1;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);

  /* allocate and initialize the workspace for rcont2. */
  workspace = alloc1dcont(nc);

  /* initialize the contingency table. */
  n = alloc1dcont(nr * nc);

  /* initialize the marginal frequencies. */
  nrowt = alloc1dcont(nr);
  ncolt = alloc1dcont(nc);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[CMC(xx[k] - 1, yy[k] - 1, nr)]++;

  /* compute the marginals. */
  for (i = 0; i < nr; i++)
    for (k = 0; k < nc; k++) {

      nrowt[i] += n[CMC(i, k, nr)];
      ncolt[k] += n[CMC(i, k, nr)];

    }/*FOR*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MUTUAL_INFORMATION:
      *observed = mc_mi(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(&nr, &nc, nrowt, ncolt, &num, fact, workspace, n);

        if (mc_mi(n, nrowt, ncolt, nr, nc, num) > *observed) {

          sequential_counter_check(*pvalue);

        }/*THEN*/

      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case PEARSON_X2:
      *observed = mc_x2(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(&nr, &nc, nrowt, ncolt, &num, fact, workspace, n);

        if (mc_x2(n, nrowt, ncolt, nr, nc, num) > *observed) {

          sequential_counter_check(*pvalue);

        }/*THEN*/

      }/*FOR*/

      break;

    case SP_MUTUAL_INFORMATION:
      *observed = mc_mi(n, nrowt, ncolt, nr, nc, num);
      *df = 0;

      for (k = 0; k < B; k++) {

        c_rcont2(&nr, &nc, nrowt, ncolt, &num, fact, workspace, n);
        *df += mc_mi(n, nrowt, ncolt, nr, nc, num);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_PEARSON_X2:
      *observed = mc_x2(n, nrowt, ncolt, nr, nc, num);
      *df = 0;

      for (k = 0; k < B; k++) {

        c_rcont2(&nr, &nc, nrowt, ncolt, &num, fact, workspace, n);
        *df += mc_x2(n, nrowt, ncolt, nr, nc, num);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case JT:
      *observed = mc_jt(n, nrowt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(&nr, &nc, nrowt, ncolt, &num, fact, workspace, n);

        if (fabs(mc_jt(n, nrowt, nr, nc, num)) >= fabs(*observed)) {

          sequential_counter_check(*pvalue);

        }/*THEN*/

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(c_jt_var(num, nrowt, nr, ncolt, nc));

      break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the p-value (for nonparametric tests) or the degrees of freedon (for
   * semiparametric tests). */
  if ((test == SP_MUTUAL_INFORMATION) || (test == SP_PEARSON_X2))
    *pvalue = pchisq(*observed, *df, FALSE, FALSE);
  else
    *pvalue /= B;


}/*C_MCARLO*/

/* conditional Monte Carlo and semiparametric discrete tests. */
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, int test, double *df) {

double *fact = NULL;
int **n = NULL, **ncolt = NULL, **nrowt = NULL, *ncond = NULL, *workspace = NULL;
int i = 0, j = 0, k = 0, enough = ceil(alpha * B) + 1;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);

  /* allocate and initialize the workspace for rcont2. */
  workspace = alloc1dcont(nc);

  /* initialize the contingency table. */
  n = alloc2dcont(nl, nr * nc);

  /* initialize the marginal frequencies. */
  nrowt = alloc2dcont(nl, nr);
  ncolt = alloc2dcont(nl, nc);
  ncond = alloc1dcont(nl);

  /* compute the joint frequency of x, y and z. */
  for (k = 0; k < num; k++)
    n[zz[k] - 1][CMC(xx[k] - 1, yy[k] - 1, nr)]++;

  /* compute the marginals. */
  for (i = 0; i < nr; i++)
    for (j = 0; j < nc; j++)
      for (k = 0; k < nl; k++) {

        nrowt[k][i] += n[k][CMC(i, j, nr)];
        ncolt[k][j] += n[k][CMC(i, j, nr)];
        ncond[k] += n[k][CMC(i, j, nr)];

      }/*FOR*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MUTUAL_INFORMATION:
      *observed = mc_cmi(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(&nr, &nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (mc_cmi(n, nrowt, ncolt, ncond, nr, nc, nl) > *observed) {

          sequential_counter_check(*pvalue);

        }/*THEN*/


      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case PEARSON_X2:
      *observed = mc_cx2(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(&nr, &nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (mc_cx2(n, nrowt, ncolt, ncond, nr, nc, nl) > *observed) {

          sequential_counter_check(*pvalue);

        }/*THEN*/

      }/*FOR*/

      break;

    case SP_MUTUAL_INFORMATION:
      *observed = mc_cmi(n, nrowt, ncolt, ncond, nr, nc, nl);
      *df = 0;

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(&nr, &nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        *df += mc_cmi(n, nrowt, ncolt, ncond, nr, nc, nl);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_PEARSON_X2:
      *observed = mc_cx2(n, nrowt, ncolt, ncond, nr, nc, nl);
      *df = 0;

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(&nr, &nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        *df += mc_cx2(n, nrowt, ncolt, ncond, nr, nc, nl);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case JT:
      *observed = mc_cjt(n, nrowt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(&nr, &nc, nrowt[k], ncolt[k], &(ncond[k]), fact, workspace, n[k]);

        if (fabs(mc_cjt(n, nrowt, ncond, nr, nc, nl)) >= fabs(*observed)) {

          sequential_counter_check(*pvalue);

        }/*THEN*/

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(mc_cvjt(nrowt, ncolt, ncond, nr, nc, nl));

      break;

  }/*SWITCH*/

  PutRNGstate();

  /* save the p-value (for nonparametric tests) or the degrees of freedon (for
   * semiparametric tests). */
  if ((test == SP_MUTUAL_INFORMATION) || (test == SP_PEARSON_X2))
    *pvalue = pchisq(*observed, *df, FALSE, FALSE);
  else
    *pvalue /= B;

}/*C_CMCARLO*/

/* compute the mutual information from the joint and marginal frequencies. */
static double mc_mi(int *n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++)
      res += MI_PART(n[CMC(i, j, nrows)], nrowt[i], ncolt[j], length);

  return res;

}/*MC_MI*/

/* compute the conditional mutual information from the joint and marginal
 * frequencies. */
static double mc_cmi(int **n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < nl; k++)
    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++)
        res += MI_PART(n[k][CMC(i, j, nr)], nrowt[k][i], ncolt[k][j], ncond[k]);

  return res;

}/*MC_CMI*/

/* compute Pearson's X^2 coefficient from the joint and marginal frequencies. */
static double mc_x2(int *n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < nrows; i++)
    for (j = 0; j < ncols; j++) {

      if (n[CMC(i, j, nrows)] != 0)
        res += (n[CMC(i, j, nrows)] - nrowt[i] * (double)ncolt[j] / length) *
               (n[CMC(i, j, nrows)] - nrowt[i] * (double)ncolt[j] / length) /
               (nrowt[i] * (double)ncolt[j] / length);

    }/*FOR*/

  return res;

}/*MC_X2*/

/* compute the Pearson's conditional X^2 coefficient from the joint and
 * marginal frequencies. */
static double mc_cx2(int **n, int **nrowt, int **ncolt, int *ncond,
    int nr, int nc, int nl) {

int i = 0, j = 0, k = 0;
double res = 0;

  for (k = 0; k < nl; k++)
    for (j = 0; j < nc; j++)
      for (i = 0; i < nr; i++) {

       if (n[k][CMC(i, j, nr)] != 0) {

          res += (n[k][CMC(i, j, nr)] - nrowt[k][i] * (double)ncolt[k][j] / ncond[k]) *
                 (n[k][CMC(i, j, nr)] - nrowt[k][i] * (double)ncolt[k][j] / ncond[k]) /
                 (nrowt[k][i] * (double)ncolt[k][j] / ncond[k]);

        }/*THEN*/

      }/*FOR*/

  return res;

}/*MC_CX2*/

/* unconditional Jonckheere-Terpstra test statistic. */
static double mc_jt(int *n, int *nrowt, int nrows, int ncols, int length) {

int i = 0, j = 0, s = 0, t = 0;
double res = 0, nrt2 = 0, wi = 0, w = 0;
double mean = c_jt_mean(length, nrowt, nrows);

  for (i = 1; i < nrows; i++) {

    nrt2 = (double)(nrowt[i]) * ((double)(nrowt[i]) + 1)/2;

    for (j = 0; j < i; j++) {

      for (s = 0, w = 0; s < ncols; s++) {

        for (t = 0, wi = 0; t < s; t++)
          wi += n[CMC(i, t, nrows)] + n[CMC(j, t, nrows)];

        w += (wi + ((double)(n[CMC(i, s, nrows)]) + (double)(n[CMC(j, s, nrows)]) + 1)/2) *
               (double)(n[CMC(i, s, nrows)]);

      }/*FOR*/

      res += w - nrt2;

    }/*FOR*/

  }/*FOR*/

  return res - mean;

}/*MC_JT*/

static double mc_cjt(int **n, int **nrowt, int *ncond, int nr, int nc, int nl) {

int k = 0;
double res = 0;

  /* sum up over the parents' configurations. */
  for (k = 0; k < nl; k++) {

    /* this one is never observed, skip. */
    if (ncond[k] == 0)
      continue;

    res += mc_jt(n[k], nrowt[k], nr, nc, ncond[k]);

  }/*FOR*/

  return res;

}/*MC_CJT*/

static double mc_cvjt(int **nrowt, int **ncolt, int *ncond, int nr, int nc,
    int nl) {

int k = 0;
double res = 0, var = 0;

  for (k = 0; k < nl; k++) {

    var = c_jt_var(ncond[k], nrowt[k], nr, ncolt[k], nc);
    if (!ISNAN(var))
      res += var;

  }/*FOR*/

  return res;

}/*MC_CVJT*/

