#include "include/rcore.h"
#include "include/globals.h"
#include "include/contingency.tables.h"
#include "include/tests.h"

static double c_jt_mean(int num, int *ni, int llx);
static double c_jt_var(int num, int *ni, int llx, int *nj, int lly);

static double c_jt_stat(int **n, int *ni, int llx, int lly);

/* unconditional Jonckheere-Terpstra asymptotic test for use in C code. */
double c_jt(int *xx, int llx, int *yy, int lly, int num) {

double stat = 0, mean = 0, var = 0, res = 0;
counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(llx, lly, TRUE);
  fill_2d_table(xx, yy, &joint, num);

  /* if there are no complete data points, or if there is just a single complete
   * observation, return independence. */
  if (joint.nobs <= 1)
    goto free_and_return;

  /* compute the test statistic, its mean and it variance. */
  stat = c_jt_stat(joint.n, joint.ni, joint.llx, joint.lly);
  mean = c_jt_mean(joint.nobs, joint.ni, joint.llx);
  var = c_jt_var(joint.nobs, joint.ni, joint.llx, joint.nj, joint.lly);

  /* standardize before returning so to make the test invariant to
   * the order of the variables; if the variance is zero return zero
   * to imply independence. */
  res = (var < MACHINE_TOL) ? 0 : (stat - mean) / sqrt(var);

free_and_return:

  Free2DTAB(joint);

  return res;

}/*C_JT*/

/* conditional Jonckheere-Terpstra asymptotic test for use in C code. */
double c_cjt(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num) {

int k = 0;
double stat = 0, var = 0, tvar = 0;
counts3d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_3d_table(llx, lly, llz);
  fill_3d_table(xx, yy, zz, &joint, num);

  /* sum up over the parents' configurations. */
  for (k = 0; k < joint.llz; k++) {

    /* this one is never observed, skip. */
    if (joint.nk[k] == 0)
      continue;

    stat += c_jt_stat(joint.n[k], joint.ni[k], joint.llx, joint.lly) -
              c_jt_mean(joint.nk[k], joint.ni[k], joint.llx);
    var = c_jt_var(joint.nk[k], joint.ni[k], joint.llx, joint.nj[k], joint.lly);
    if (!ISNAN(var) && (var > MACHINE_TOL))
      tvar += var;

  }/*FOR*/

  Free3DTAB(joint);

  /* guard against divide-by-zero errors and standardize the result. */
  return (tvar < MACHINE_TOL) ? 0 : stat / sqrt(tvar);

}/*C_CJT*/

/* asymptotic expected value of the Jonkheere-Terpstra test statistic. */
static double c_jt_mean(int num, int *ni, int llx) {

double res = (double)(num * num);

  for (int i = 0; i < llx; i++)
    res -= ni[i] * ni[i];

  return res / 4;

}/*C_JT_MEAN*/

/* asymptotic variance of the Jonkheere-Terpstra test statistic. */
static double c_jt_var(int num, int *ni, int llx, int *nj, int lly) {

int i = 0, j = 0;
double U1 = 0, U2a = 0, U2b = 0, U2 = 0, U3a = 0, U3b = 0, U3 = 0;
double res = 0, t1 = 0, t2 = 0, t3 = 0;

  /* numerator, first term. */
  U1 = (double)(num) * (double)(num - 1) * (double)(2 * num + 5);
  for (i = 0; i < llx; i++)
    U1 -= (double)(ni[i]) * (double)(ni[i] -1) * (double)(2 * ni[i] + 5);
  for (j = 0; j < lly; j++)
    U1 -= (double)(nj[j]) * (double)(nj[j] -1) * (double)(2 * nj[j] + 5);

  /* numerator, second term. */
  for (i = 0; i < llx; i++)
    U2a += (double)(ni[i]) * (double)(ni[i]  - 1) * (double)(ni[i] - 2);
  for (j = 0; j < lly; j++)
    U2b += (double)(nj[j]) * (double)(nj[j]  - 1) * (double)(nj[j] - 2);
  U2 = U2a * U2b;

  /* numerator, third term. */
  for (i = 0; i < llx; i++)
    U3a += (double)(ni[i]) * (double)(ni[i]  - 1);
  for (j = 0; j < lly; j++)
    U3b += (double)(nj[j]) * (double)(nj[j]  - 1);
  U3 = U3a * U3b;

  /* terms in the denominators. */
  t1 = 72;
  t2 = 36 * (double)(num) * (double)(num - 1) * (double)(num - 2);
  t3 = 8 * (double)(num) * (double)(num - 1);

  res = U1/t1 + U2/t2 + U3/t3;

  return res;

}/*C_JT_VAR*/

/* unconditional Jonckheere-Terpstra test statistic. */
static double c_jt_stat(int **n, int *ni, int llx, int lly) {

int i = 0, j = 0, s = 0, t = 0;
double res = 0, ni2 = 0, wi = 0, w = 0;

  for (i = 1; i < llx; i++) {

    ni2 = (double)(ni[i]) * ((double)(ni[i]) + 1)/2;

    for (j = 0; j < i; j++) {

      for (s = 0, w = 0; s < lly; s++) {

        for (t = 0, wi = 0; t < s; t++)
          wi += n[i][t] + n[j][t];

        w += (wi + ((double)(n[i][s]) + (double)(n[j][s]) + 1)/2) * (double)(n[i][s]);

      }/*FOR*/

      res += w - ni2;

    }/*FOR*/

  }/*FOR*/

  return res;

}/*C_JT_STAT*/

/* compute the centered Jonckheere-Terpstra statistic. */
double jt_centered_kernel(counts2d table) {

  return c_jt_stat(table.n, table.ni, table.llx, table.lly) -
           c_jt_mean(table.nobs, table.ni, table.llx);

}/*JT_CENTERED_KERNEL*/

/* compute the conditional centered Jonckheere-Terpstra statistic. */
double cjt_centered_kernel(counts3d table) {

long double res = 0;

  /* sum up over the parents' configurations. */
  for (int k = 0; k < table.llz; k++) {

    /* this one is never observed, skip. */
    if (table.nk[k] == 0)
       continue;

    res += c_jt_stat(table.n[k], table.ni[k], table.llx, table.lly) -
             c_jt_mean(table.nk[k], table.ni[k], table.llx);

  }/*FOR*/

  return (double)res;

}/*CJT_CENTERED_KERNEL*/

/* compute the asymptotic variance of the Jonckheere-Terpstra statistic. */
double jt_var_kernel(counts2d table) {

  return c_jt_var(table.nobs, table.ni, table.llx, table.nj, table.lly);

}/*JT_VAR_KERNEL*/

/* compute the asymptotic variance of the conditional Jonckheere-Terpstra
 * statistic. */
double cjt_var_kernel(counts3d table) {

long double res = 0;
double var = 0;

  for (int k = 0; k < table.llz; k++) {

    var = c_jt_var(table.nk[k], table.ni[k], table.llx, table.nj[k], table.lly);
    if (!ISNAN(var))
      res += var;

  }/*FOR*/

  return (double)res;

}/*CJT_VAR_KERNEL*/

