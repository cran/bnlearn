#include "include/rcore.h"
#include "include/sampling.h"
#include "include/tests.h"
#include "include/matrix.h"

/* initialize the table of log-factorials. */
#define allocfact(n) \
  fact = Calloc1D(n + 1, sizeof(double)); \
  fact[0] = 0.; \
  for(k = 1; k <= n; k++) \
    fact[k] = lgammafn((double) (k + 1));

static double mc_jt(int **n, int *nrowt, int nrows, int ncols, int length);
static double mc_cjt(int ***n, int **nrowt, int *ncond, int nr, int nc, int nl);
static double mc_cvjt(int **nrowt, int **ncolt, int *ncond, int nr, int nc,
    int nl);

/* unconditional Monte Carlo and semiparametric discrete tests. */
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, test_e test, double *df) {

double *fact = NULL;
int **n = NULL, *ncolt = NULL, *nrowt = NULL, *workspace = NULL;
int k = 0, enough = ceil(alpha * B) + 1, constx = TRUE, consty = TRUE;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);
  /* allocate and initialize the workspace for rcont2. */
  workspace = Calloc1D(nc, sizeof(int));
  /* initialize the contingency table and the marginal frequencies. */
  fill_2d_table(xx, yy, &n, &nrowt, &ncolt, nr, nc, num);

  /* if at least one of the two variables is constant, they are independent. */
  for (k = 0; k < nr; k++)
    constx = constx && ((nrowt[k] == 0) || (nrowt[k] == num));
  for (k = 0; k < nc; k++)
    consty = consty && ((ncolt[k] == 0) || (ncolt[k] == num));

  if (constx || consty) {

    *observed = 0;
    *pvalue = 1;

    Free2D(n, nr);
    Free1D(nrowt);
    Free1D(ncolt);
    Free1D(fact);
    Free1D(workspace);

    return;

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MC_MI:
    case SMC_MI:
      *observed = mi_kernel(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (mi_kernel(n, nrowt, ncolt, nr, nc, num) > *observed) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/

      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case MC_X2:
    case SMC_X2:
      *observed = x2_kernel(n, nrowt, ncolt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (x2_kernel(n, nrowt, ncolt, nr, nc, num) > *observed) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/

      }/*FOR*/

      break;

    case SP_MI:
      *observed = mi_kernel(n, nrowt, ncolt, nr, nc, num);
      *df = 0;

      for (k = 0; k < B; k++) {

        c_rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);
        *df += mi_kernel(n, nrowt, ncolt, nr, nc, num);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_X2:
      *observed = x2_kernel(n, nrowt, ncolt, nr, nc, num);
      *df = 0;

      for (k = 0; k < B; k++) {

        c_rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);
        *df += x2_kernel(n, nrowt, ncolt, nr, nc, num);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case MC_JT:
    case SMC_JT:
      *observed = mc_jt(n, nrowt, nr, nc, num);

      for (k = 0; k < B; k++) {

        c_rcont2(nr, nc, nrowt, ncolt, num, fact, workspace, n);

        if (fabs(mc_jt(n, nrowt, nr, nc, num)) >= fabs(*observed)) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(c_jt_var(num, nrowt, nr, ncolt, nc));

      break;

    default:
      error("unknown permutation test statistic.");

  }/*SWITCH*/

  PutRNGstate();

  /* save the p-value (for nonparametric tests) or the degrees of freedon (for
   * semiparametric tests). */
  if ((test == SP_MI) || (test == SP_X2))
    *pvalue = pchisq(*observed, *df, FALSE, FALSE);
  else
    *pvalue /= B;

  Free2D(n, nr);
  Free1D(nrowt);
  Free1D(ncolt);
  Free1D(workspace);
  Free1D(fact);

}/*C_MCARLO*/

/* conditional Monte Carlo and semiparametric discrete tests. */
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, test_e test,
    double *df) {

double *fact = NULL;
int ***n = NULL, **ncolt = NULL, **nrowt = NULL, *ncond = NULL, *workspace = NULL;
int j = 0, k = 0, enough = ceil(alpha * B) + 1, constx = TRUE, consty = TRUE;

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);
  /* allocate and initialize the workspace for rcont2. */
  workspace = Calloc1D(nc, sizeof(int));
  /* initialize the contingency table and the marginal frequencies. */
  fill_3d_table(xx, yy, zz, &n, &nrowt, &ncolt, &ncond, nr, nc, nl, num);

  /* if at least one of the two variables is constant, they are independent. */
  for (k = 0; k < nl; k++)
    for (j = 0; j < nr; j++)
      constx = constx && ((nrowt[k][j] == 0) || (nrowt[k][j] == ncond[k]));
  for (k = 0; k < nl; k++)
    for (j = 0; j < nc; j++)
      consty = consty && ((ncolt[k][j] == 0) || (ncolt[k][j] == ncond[k]));

  if (constx || consty) {

    *observed = 0;
    *pvalue = 1;

    Free3D(n, nl, nr);
    Free2D(nrowt, nl);
    Free2D(ncolt, nl);
    Free1D(ncond);
    Free1D(fact);
    Free1D(workspace);

    return;

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MC_MI:
    case SMC_MI:
      *observed = cmi_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(nr, nc, nrowt[k], ncolt[k], ncond[k], fact, workspace, n[k]);

        if (cmi_kernel(n, nrowt, ncolt, ncond, nr, nc, nl) > *observed) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/


      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case MC_X2:
    case SMC_X2:
      *observed = cx2_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(nr, nc, nrowt[k], ncolt[k], ncond[k], fact, workspace, n[k]);

        if (cx2_kernel(n, nrowt, ncolt, ncond, nr, nc, nl) > *observed) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/

      }/*FOR*/

      break;

    case SP_MI:
      *observed = cmi_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);
      *df = 0;

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(nr, nc, nrowt[k], ncolt[k], ncond[k], fact, workspace, n[k]);

        *df += cmi_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_X2:
      *observed = cx2_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);
      *df = 0;

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(nr, nc, nrowt[k], ncolt[k], ncond[k], fact, workspace, n[k]);

        *df += cx2_kernel(n, nrowt, ncolt, ncond, nr, nc, nl);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case MC_JT:
    case SMC_JT:
      *observed = mc_cjt(n, nrowt, ncond, nr, nc, nl);

      for (j = 0; j < B; j++) {

        for (k = 0; k < nl; k++)
          c_rcont2(nr, nc, nrowt[k], ncolt[k], ncond[k], fact, workspace, n[k]);

        if (fabs(mc_cjt(n, nrowt, ncond, nr, nc, nl)) >= fabs(*observed)) {

          SEQUENTIAL_COUNTER_CHECK(*pvalue);

        }/*THEN*/

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(mc_cvjt(nrowt, ncolt, ncond, nr, nc, nl));

      break;

    default:
      error("unknown permutation test statistic.");

  }/*SWITCH*/

  PutRNGstate();

  /* save the p-value (for nonparametric tests) or the degrees of freedon (for
   * semiparametric tests). */
  if ((test == SP_MI) || (test == SP_X2))
    *pvalue = pchisq(*observed, *df, FALSE, FALSE);
  else
    *pvalue /= B;

  Free3D(n, nl, nr);
  Free2D(nrowt, nl);
  Free2D(ncolt, nl);
  Free1D(ncond);
  Free1D(workspace);
  Free1D(fact);

}/*C_CMCARLO*/

/* unconditional Jonckheere-Terpstra test statistic. */
static double mc_jt(int **n, int *nrowt, int nrows, int ncols, int length) {

int i = 0, j = 0, s = 0, t = 0;
double res = 0, nrt2 = 0, wi = 0, w = 0;
double mean = c_jt_mean(length, nrowt, nrows);

  for (i = 1; i < nrows; i++) {

    nrt2 = (double)(nrowt[i]) * ((double)(nrowt[i]) + 1)/2;

    for (j = 0; j < i; j++) {

      for (s = 0, w = 0; s < ncols; s++) {

        for (t = 0, wi = 0; t < s; t++)
          wi += n[i][t] + n[j][t];

        w += (wi + ((double)(n[i][s]) + (double)(n[j][s]) + 1)/2) *
               (double)(n[i][s]);

      }/*FOR*/

      res += w - nrt2;

    }/*FOR*/

  }/*FOR*/

  return res - mean;

}/*MC_JT*/

static double mc_cjt(int ***n, int **nrowt, int *ncond, int nr, int nc, int nl) {

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

