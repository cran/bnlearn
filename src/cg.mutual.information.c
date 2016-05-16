#include "include/rcore.h"
#include "include/scores.h"
#include "include/sets.h"
#include "include/tests.h"

/* unconditional mutual information, to be used in C code. */
double c_micg(double *yy, double ym, double ysd, int *xx, int llx, int num) {

int i = 0, *ni = NULL;
double lognum = 0, logden = 0, *mu = NULL, *sd = NULL;

  /* compute the denominator (model under the null). */
  for (i = 0; i < num; i++)
    logden += dnorm(yy[i], ym, ysd, TRUE);

  mu = Calloc1D(llx, sizeof(double));
  sd = Calloc1D(llx, sizeof(double));
  ni = Calloc1D(llx, sizeof(int));

  /* compute the conditional means and variances. */
  for (i = 0; i < num; i++) {

    mu[xx[i] - 1] += yy[i];
    ni[xx[i] - 1]++;

  }/*FOR*/
  for (i = 0; i < llx; i++)
    mu[i] /= ni[i];

  for (i = 0; i < num; i++)
    sd[xx[i] - 1] += (yy[i] - mu[xx[i] - 1]) * (yy[i] - mu[xx[i] - 1]);
  for (i = 0; i < llx; i++)
    sd[i] = sqrt(sd[i] / (ni[i] - 1));

  /* compute the numerator (model under the alternative). */
  for (i = 0; i < num; i++)
    lognum += dnorm(yy[i], mu[xx[i] - 1], sd[xx[i] - 1], TRUE);

  Free1D(mu);
  Free1D(sd);
  Free1D(ni);

  return (lognum - logden) / num;

}/*C_MICG*/

/* conditional mutual information, to be used in C code. */
double c_cmicg(double *yy, double **xx, int nx, int **zz, int nz, int *z0,
    int nz0, int *nlvls, int num) {

double logden = 0, lognum = 0;
int *z1 = NULL, nz1 = 0;

  if (!zz) {

    /* compute the denominator (model under the null). */
    logden = c_fast_ccgloglik(yy, xx + 1, nx - 1, num, z0, nz0);
    /* compute the numerator (model under the alternative). */
    lognum = c_fast_ccgloglik(yy, xx, nx, num, z0, nz0);

  }/*THEN*/
  else {

    /* compute the denominator (model under the null). */
    logden = c_fast_ccgloglik(yy, xx, nx, num, z0, nz0);
    /* compute the numerator (model under the alternative). */
    z1 = Calloc1D(num, sizeof(int));
    c_fast_config(zz, num, nz, nlvls, z1, &nz1, 1);
    lognum = c_fast_ccgloglik(yy, xx, nx, num, z1, nz1);
    Free1D(z1);

  }/*ELSE*/

  /* if the null model is singular, the alternative model is even more singular
   * so it should always be rejected. */
  return (R_FINITE(logden) && R_FINITE(lognum)) ? (lognum - logden) / num : 0;

}/*C_CMIGG*/

/* conditional mutual information between two discrete variables, conditional
 * on both discrete and continuous variables. */
double c_cmicg_unroll(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    double **gp, int ngp, double *df, int num) {

int i = 0, llz2 = 0, tlvls[2] = {0, 0};
int *tt[2] = {NULL, NULL}, *zz2 = NULL;
double logden = 0, lognum = 0;

  if (!zz) {

    /* no discrete conditioning variables, reuse yy.n the denominator. */
    zz2 = yy;
    llz2 = lly;

    /* remainder term: the mutual information test of xx and yy. */
    lognum = c_chisqtest(xx, llx, yy, lly, num, df, MI);

  }/*THEN*/
  else {

    /* combine zz and yy to get the configurations for the denominator. */
    zz2 = Calloc1D(num, sizeof(int));
    tt[0] = yy;
    tt[1] = zz;
    tlvls[0] = lly;
    tlvls[1] = llz;
    c_fast_config(tt, num, 2, tlvls, zz2, &llz2, 1);

    /* remainder term: the mutual information test of xx and yy given zz. */
    lognum = c_cchisqtest(xx, llx, yy, lly, zz, llz, num, df, MI);

  }/*ELSE*/

  /* iterate over the continuous variables using the chain rule; using only the
   * first variable xx for the numerator and both xx and yy for the denomiator. */
  for (i = 0; i < ngp; i++)
    lognum += c_cmicg(gp[i], gp + i + 1, ngp - i - 1, &xx , 1, zz, llz, &llx, num);

  for (i = 0; i < ngp; i++)
    logden += c_cmicg(gp[i], gp + i + 1, ngp - i - 1, &xx , 1, zz2, llz2, &llx, num);

  if (zz)
    Free1D(zz2);

  /* set the degrees of freedom. */
  if (df)
    *df += llz * (ngp * (ngp + 3) / 2) * (lly - 1);

  /* if the null model is singular, the alternative model is even more singular
   * so it should always be rejected. */
  return (R_FINITE(logden) && R_FINITE(lognum)) ? (lognum - logden) / num : 0;

}/*C_CMICG_UNROLL*/

