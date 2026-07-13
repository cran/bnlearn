#include "../include/rcore.h"
#include "../include/globals.h"
#include "./xnegbin.h"

/* density of an real-valued negative binomial random variable. */
double dxnegbin(double x, double p, double fails, bool give_log) {

double logf = 0;

  if (ISNAN(x) || ISNAN(p) || ISNAN(fails))
    return NA_REAL;

  logf = fails * log(1 - p);

  if (x > 0)
    logf += lgammafn(fails + x) - lgammafn(fails) - lgammafn(x + 1) +
              x * log(p);

  if (give_log)
    return logf;
  else
    return exp(logf);

}/*DXNEGBIN*/

/* mean of a real-valued negative binomial random variable: the NB2-parameterised
 * mean fails * p / (1 - p). */
double negbin_mean(double p, double fails) {

  return fails * p / (1 - p);

}/*NEGBIN_MEAN*/

/* random sampling from a real-valued negative binomial random variable. */
double rxnegbin(double p, double fails) {

  if (ISNAN(p) || ISNAN(fails))
    return NA_REAL;
  if (fails < MACHINE_TOL)
    return 0;
  if (p < MACHINE_TOL)
    return 0;
  if (p > 1 - MACHINE_TOL)
    return R_PosInf;

  /* Use the gamma-Poisson mixture representation, as R does. */
  return rpois(rgamma(fails, p / (1 - p)));

}/*RXNEGBIN*/
