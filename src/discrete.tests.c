#include "include/rcore.h"
#include "include/contingency.tables.h"
#include "include/tests.h"

#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))

/* unconditional parametric asymptotic tests for categorical data. */
double c_chisqtest(int *xx, int llx, int *yy, int lly, int num, double *df,
    test_e test, bool scale) {

double res = 0;
counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(llx, lly, TRUE);
  fill_2d_table(xx, yy, &joint, num);

  /* compute the degrees of freedom. */
  if (df)
    *df = discrete_df(test, joint.ni, joint.llx, joint.nj, joint.lly);

  /* if there are no complete data points, return independence. */
  if (joint.nobs == 0)
    goto free_and_return;

  /* if there are less than 5 observations per cell on average, assume the
   * test does not have enough power and return independence. */
  if ((test == MI_ADF) || (test == X2_ADF))
    if (joint.nobs < 5 * joint.llx * joint.lly)
      goto free_and_return;

  /* compute the mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = mi_kernel(joint) / joint.nobs;
  else if ((test == X2) || (test == X2_ADF))
    res = x2_kernel(joint);

  /* rescale to match the G^2 test. */
  if (scale)
    res *= 2 * joint.nobs;

free_and_return:

  Free2DTAB(joint);

  return res;

}/*C_CHISQTEST*/

/* conditional mutual information, to be used in C code. */
double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, bool scale) {

double res = 0;
counts3d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_3d_table(llx, lly, llz);
  fill_3d_table(xx, yy, zz, &joint, num);

  /* compute the degrees of freedom. */
  if (df)
    *df = discrete_cdf(test, joint.ni, joint.llx, joint.nj, joint.lly, joint.llz);

  /* if there are no complete data points, return independence. */
  if (joint.nobs == 0)
    goto free_and_return;

  /* if there are less than 5 observations per cell on average, assume the
   * test does not have enough power and return independence. */
  if ((test == MI_ADF) || (test == X2_ADF))
    if (joint.nobs < 5 * joint.llx * joint.lly * joint.llz)
      goto free_and_return;

  /* compute the conditional mutual information or Pearson's X^2. */
  if ((test == MI) || (test == MI_ADF))
    res = cmi_kernel(joint) / joint.nobs;
  else if ((test == X2) || (test == X2_ADF))
    res = cx2_kernel(joint);

  /* rescale to match the G^2 test. */
  if (scale)
    res *= 2 * joint.nobs;

free_and_return:

  Free3DTAB(joint);

  return res;

}/*C_CCHISQTEST*/

/* compute the mutual information. */
double mi_kernel(counts2d table) {

long double res = 0;

  for (int i = 0; i < table.llx; i++)
    for (int j = 0; j < table.lly; j++)
      res += MI_PART(table.n[i][j], table.ni[i], table.nj[j], table.nobs);

  return (double)res;

}/*MI_KERNEL*/

/* midified mutual information for Hartemink's discretization. */
double mi_kernel_collapsed(counts2d table, int k) {

long double res = 0;

  for (int i = 0; i < table.llx; i++) {

    if ((i == k) || (i == k + 1))
      continue;

    for (int j = 0; j < table.lly; j++)
      res += MI_PART(table.n[i][j], table.ni[i], table.nj[j], table.nobs);

  }/*FOR*/

  for (int j = 0; j < table.lly; j++)
    res += MI_PART((table.n[k][j] + table.n[k + 1][j]),
             table.ni[k] + table.ni[k + 1], table.nj[j], table.nobs);

  return (double)res;

}/*MI_KERNEL_COLLAPSED*/

/* compute Pearson's X^2 coefficient. */
double x2_kernel(counts2d table) {

long double res = 0, expected = 0;

  for (int i = 0; i < table.llx; i++)
    for (int j = 0; j < table.lly; j++) {

      expected = table.ni[i] * (double)table.nj[j] / table.nobs;

      if (expected != 0)
        res += (table.n[i][j] - expected) *
               (table.n[i][j] - expected) / expected;

    }/*FOR*/

  return (double)res;

}/*X2_KERNEL*/

/* compute the conditional mutual information. */
double cmi_kernel(counts3d table) {

long double res = 0;

  for (int k = 0; k < table.llz; k++)
    for (int i = 0; i < table.llx; i++)
      for (int j = 0; j < table.lly; j++)
        res += MI_PART(table.n[k][i][j], table.ni[k][i], table.nj[k][j],
                 table.nk[k]);

  return (double)res;

}/*CMI_KERNEL*/

/* compute the Pearson's conditional X^2 coefficient. */
double cx2_kernel(counts3d table) {

long double expected = 0, res = 0;

  for (int k = 0; k < table.llz; k++) {

    if (table.nk[k] == 0)
      continue;

    for (int i = 0; i < table.llx; i++)
      for (int j = 0; j < table.lly; j++) {

       expected = table.ni[k][i] * (double)table.nj[k][j] / table.nk[k];

       if (expected != 0)
          res += (table.n[k][i][j] - expected) *
                 (table.n[k][i][j] - expected) / expected;

      }/*FOR*/

  }/*FOR*/

  return (double)res;

}/*CX2_KERNEL*/

