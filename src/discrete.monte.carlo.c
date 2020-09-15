#include "include/rcore.h"
#include "include/sampling.h"
#include "include/contingency.tables.h"
#include "include/tests.h"
#include "include/matrix.h"

/* initialize the table of log-factorials. */
#define allocfact(n) \
  fact = Calloc1D(n + 1, sizeof(double)); \
  fact[0] = 0.; \
  for(k = 1; k <= n; k++) \
    fact[k] = lgammafn((double) (k + 1));

/* unconditional Monte Carlo and semiparametric discrete tests. */
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, test_e test, double *df) {

double *fact = NULL;
int *workspace = NULL;
int k = 0, enough = ceil(alpha * B) + 1, constx = TRUE, consty = TRUE;
counts2d joint = { 0 };

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);
  /* allocate and initialize the workspace for rcont2. */
  workspace = Calloc1D(nc, sizeof(int));
  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(nr, nc, TRUE);
  fill_2d_table(xx, yy, &joint, num);

  /* if at least one of the two variables is constant, or if there are no
   * complete observations, the variables are taken to be independent. */
  for (k = 0; k < joint.llx; k++)
    constx = constx && ((joint.ni[k] == 0) || (joint.ni[k] == joint.nobs));
  for (k = 0; k < joint.lly; k++)
    consty = consty && ((joint.nj[k] == 0) || (joint.nj[k] == joint.nobs));

  if (constx || consty) {

    *observed = 0;
    *pvalue = 1;

    goto free_and_return;

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MC_MI:
    case SMC_MI:
      *observed = mi_kernel(joint);

      for (k = 0, *pvalue = 0; k < B; k++) {

        rcounts2d(joint, fact, workspace);
        if (mi_kernel(joint) >= *observed)
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case MC_X2:
    case SMC_X2:
      *observed = x2_kernel(joint);

      for (k = 0, *pvalue = 0; k < B; k++) {

        rcounts2d(joint, fact, workspace);
        if (x2_kernel(joint) >= *observed)
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      break;

    case SP_MI:
      *observed = mi_kernel(joint);
      *df = 0;

      for (k = 0, *pvalue = 0; k < B; k++) {

        rcounts2d(joint, fact, workspace);
        *df += mi_kernel(joint);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_X2:
      *observed = x2_kernel(joint);
      *df = 0;

      for (k = 0, *pvalue = 0; k < B; k++) {

        rcounts2d(joint, fact, workspace);
        *df += x2_kernel(joint);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case MC_JT:
    case SMC_JT:
      *observed = jt_centered_kernel(joint);

      for (k = 0, *pvalue = 0; k < B; k++) {

        rcounts2d(joint, fact, workspace);

        if (fabs(jt_centered_kernel(joint)) >= fabs(*observed))
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(jt_var_kernel(joint));

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

free_and_return:

  Free1D(workspace);
  Free1D(fact);
  Free2DTAB(joint);

}/*C_MCARLO*/

/* conditional Monte Carlo and semiparametric discrete tests. */
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, test_e test,
    double *df) {

double *fact = NULL;
int *workspace = NULL;
int j = 0, k = 0, enough = ceil(alpha * B) + 1, constx = TRUE, consty = TRUE;
counts3d joint = { 0 };

  /* allocate and compute the factorials needed by rcont2. */
  allocfact(num);
  /* allocate and initialize the workspace for rcont2. */
  workspace = Calloc1D(nc, sizeof(int));
  /* initialize the contingency table and the marginal frequencies. */
  joint = new_3d_table(nr, nc, nl);
  fill_3d_table(xx, yy, zz, &joint, num);

  /* if at least one of the two variables is constant, or if there are no
   * complete observations, the variables are taken to be independent. */
  for (k = 0; k < joint.llz; k++)
    for (j = 0; j < joint.llx; j++)
      constx = constx && ((joint.ni[k][j] == 0) || (joint.ni[k][j] == joint.nk[k]));
  for (k = 0; k < joint.llz; k++)
    for (j = 0; j < joint.lly; j++)
      consty = consty && ((joint.nj[k][j] == 0) || (joint.nj[k][j] == joint.nk[k]));

  if (constx || consty) {

    *observed = 0;
    *pvalue = 1;

    goto free_and_return;

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* pick up the observed value of the test statistic, then generate a set of
     random contingency tables (given row and column totals) and check how many
     tests are greater than the original one.*/
  switch(test) {

    case MC_MI:
    case SMC_MI:
      *observed = cmi_kernel(joint);

      for (j = 0, *pvalue = 0; j < B; j++) {

        rcounts3d(joint, fact, workspace);
        if (cmi_kernel(joint) >= *observed)
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      *observed = 2 * (*observed);

      break;

    case MC_X2:
    case SMC_X2:
      *observed = cx2_kernel(joint);

      for (j = 0, *pvalue = 0; j < B; j++) {

        rcounts3d(joint, fact, workspace);
        if (cx2_kernel(joint) >= *observed)
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      break;

    case SP_MI:
      *observed = cmi_kernel(joint);
      *df = 0;

      for (j = 0, *pvalue = 0; j < B; j++) {

        rcounts3d(joint, fact, workspace);
        *df += cmi_kernel(joint);

      }/*FOR*/

      *observed = 2 * (*observed);

      /* estimate the degrees of freedom as the expectation under the null. */
      *df = (*df) * 2 / B;

      break;

    case SP_X2:
      *observed = cx2_kernel(joint);
      *df = 0;

      for (j = 0, *pvalue = 0; j < B; j++) {

        rcounts3d(joint, fact, workspace);
        *df += cx2_kernel(joint);

      }/*FOR*/

      /* estimate the degrees of freedom as the expectation under the null. */
      *df /= B;

      break;

    case MC_JT:
    case SMC_JT:
      *observed = cjt_centered_kernel(joint);

      for (j = 0, *pvalue = 0; j < B; j++) {

        rcounts3d(joint, fact, workspace);

        if (fabs(cjt_centered_kernel(joint)) >= fabs(*observed))
          SEQUENTIAL_COUNTER_CHECK(*pvalue);

      }/*FOR*/

      /* standardize to match the parametric test. */
      *observed /= sqrt(cjt_var_kernel(joint));

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

free_and_return:

  Free3DTAB(joint);
  Free1D(workspace);
  Free1D(fact);

}/*C_CMCARLO*/

