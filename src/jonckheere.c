#include "include/rcore.h"
#include "include/globals.h"
#include "include/allocations.h"
#include "include/tests.h"

static double c_jt_stat(int **n, int *ni, int llx, int lly);

/* unconditional Jonckheere-Terpstra asymptotic test for use in C code. */
double c_jt(int *xx, int llx, int *yy, int lly, int num) {

int i = 0, j = 0, k = 0;
int  **n = NULL, *ni = NULL, *nj = NULL;
double stat = 0, mean = 0, var = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(llx, lly);
  ni = alloc1dcont(llx);
  nj = alloc1dcont(lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* compute the test statistic, its mean and it variance. */
  stat = c_jt_stat(n, ni, llx, lly);
  mean = c_jt_mean(num, ni, llx);
  var = c_jt_var(num, ni, llx, nj, lly);

  /* standardize before returning so to make the test invariant to
   * the order of the variables; if the variance is zero return zero
   * to imply independence. */
  return (var < MACHINE_TOL) ? 0 : (stat - mean) / sqrt(var);

}/*C_JT*/

/* conditional Jonckheere-Terpstra asymptotic test for use in C code. */
double c_cjt(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num) {

int i = 0, j = 0, k = 0;
int  ***n = NULL, **nrowt = NULL, **ncolt = NULL, *ncond = NULL;
double stat = 0, var = 0, tvar = 0;

  /* initialize the contingency table. */
  n = alloc3dcont(llz, llx, lly);

  /* initialize the marginal frequencies. */
  nrowt = alloc2dcont(llz, llx);
  ncolt = alloc2dcont(llz, lly);
  ncond = alloc1dcont(llz);

  /* compute the joint frequencies. */
  for (k = 0; k < num; k++)
    n[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        nrowt[k][i] += n[k][i][j];
        ncolt[k][j] += n[k][i][j];
        ncond[k] += n[k][i][j];

      }/*FOR*/

  /* sum up over the parents' configurations. */
  for (k = 0; k < llz; k++) {

    /* this one is never observed, skip. */
    if (ncond[k] == 0) continue;

    stat += c_jt_stat(n[k], nrowt[k], llx, lly) -
              c_jt_mean(ncond[k], nrowt[k], llx);
    var = c_jt_var(ncond[k], nrowt[k], llx, ncolt[k], lly);
    if (!ISNAN(var) && (var > MACHINE_TOL))
      tvar += var;

  }/*FOR*/

  /* guard against divide-by-zero errors and standardize the result. */
  return (tvar < MACHINE_TOL) ? 0 : stat / sqrt(tvar);

}/*C_CJT*/

/* asymptotic expected value of the Jonkheere-Terpstra test statistic. */
double c_jt_mean(int num, int *ni, int llx) {

int i = 0;
double res = (double)(num * num);

  for (i = 0; i < llx; i++)
    res -= ni[i] * ni[i];

  return res / 4;

}/*C_JT_MEAN*/

/* asymptotic variance of the Jonkheere-Terpstra test statistic. */
double c_jt_var(int num, int *ni, int llx, int *nj, int lly) {

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
