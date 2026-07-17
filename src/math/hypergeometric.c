#include "../include/rcore.h"
#include "./hypergeometric.h"

/* below this z the 1F1 series is summed in ordinary (linear) space with a running
 * product rather than term-by-term in log space: the terms are all positive for
 * the (a, b, z > 0) call sites so there is no cancellation, long double keeps the
 * partial sums (which grow like e^z) finite with a wide margin, and the result is
 * exact while avoiding a logspace_add and several log calls per term. */
#define LIN_1F1_MAX 700.0

/* above this z, and provided z >> b (z > K b, keeping clear of the coalescing
 * b ~ z region where the expansion is invalid), 1F1 is taken from its large-
 * argument asymptotic directly rather than by summing the O(z)-term series: for
 * a = 1 and a = 2 the algebraic part of the expansion terminates, so the error
 * is the exponentially small recessive term, negligible for z past this bound. */
#define ASYMP_1F1_MIN 50.0
#define ASYMP_1F1_K   8.0

/* confluent hypergeometric function 1F1(a; b; z), on the log scale. */
double log1F1(double a, double b, double z) {

double a0 = a, b0 = b;
long double fac = 0, temp = 0, series = 0;

  /* large z with z >> b: use the O(1) large-argument asymptotic directly. */
  if (z > ASYMP_1F1_MIN && z > ASYMP_1F1_K * b)
    return lgammafn(b0) - lgammafn(a0) + z + (a0 - b0) * log(z) +
             log1p((b0 - a0) * (1 - a0) / z);

  /* moderate z: exact linear-space running-product summation. */
  if (z >= 0 && z < LIN_1F1_MAX) {

    long double t = 1, sum = 1;

    for (int n = 1; n < MAXITER_1F1 + 1; n++) {

      t *= (long double)(a + n - 1) * z / ((long double)(b + n - 1) * n);
      sum += t;

      if (t <= sum * TOL_1F1)
        return (double) logl(sum);

    }/*FOR*/

  }/*THEN*/

  /* large z: log-space series with a large-argument asymptotic fallback. */
  for (int n = 1; n < MAXITER_1F1 + 1; n++) {

    fac += log(a) - log(b) + log(z) - log(n);

    if (ISNAN(fac))
      fac = 0;

    series = logspace_add(temp, fac);

    /* the increment series - temp is already a (non-negative) log-scale
     * quantity, so it is compared against TOL_1F1 directly. */
    if (!R_FINITE(series) || (fabsl(series - temp) < TOL_1F1))
      return series;

    temp = series;
    a++;
    b++;

  }/*FOR*/

  /* the series peaks around n = z and so does not converge within MAXITER_1F1
   * iterations for large z; fall back to the large-argument asymptotic
   * 1F1(a; b; z) ~ Gamma(b) / Gamma(a) * exp(z) * z^(a - b) * (1 + (b - a)(1 - a) / z),
   * whose single correction term makes it exact for a = 1 and a = 2 (the only
   * values used at any call site). */
  return lgammafn(b0) - lgammafn(a0) + z + (a0 - b0) * log(z) +
           log1p((b0 - a0) * (1 - a0) / z);

}/*LOG1F1*/

/* derivative of log 1F1(1; b; z) with respect to b, computed in log-space so it
 * stays finite for large z (where 1F1 itself overflows). */
double dlog1F1_db(double b, double z) {

double b0 = b, logz = log(z), psi = 0;
long double logt = 0, logden = 0, lognum = R_NegInf, newden = 0;

  /* large z with z >> b: use the O(1) digamma asymptotic directly (see log1F1). */
  if (z > ASYMP_1F1_MIN && z > ASYMP_1F1_K * b)
    return digamma(b0) - logz;

  /* moderate z: exact linear-space running-product summation (see log1F1). */
  if (z >= 0 && z < LIN_1F1_MAX) {

    long double t = 1, den = 1, num = 0, ps = 0;

    for (int n = 1; n < MAXITER_1F1 + 1; n++) {

      ps += 1 / (long double)(b + n - 1);   /* psi_n = sum_{k < n} 1 / (b + k). */
      t *= z / (long double)(b + n - 1);     /* t_n = z^n / (b)_n. */
      num += ps * t;
      den += t;

      if (t <= den * TOL_1F1)
        return -(double)(num / den);

    }/*FOR*/

  }/*THEN*/

  /* large z: log-space series with a digamma asymptotic fallback. */
  for (int n = 1; n < MAXITER_1F1 + 1; n++) {

    /* advance log(t_n) = log(z^n / (b)_n) and psi_n = sum_{k < n} 1 / (b + k);
     * d/db log 1F1 = -sum_{n >= 1} psi_n t_n / sum_{n >= 0} t_n. */
    logt += logz - log(b);
    psi += 1 / b;

    lognum = logspace_add(lognum, log(psi) + logt);
    newden = logspace_add(logden, logt);

    if (!R_FINITE(newden) || (fabsl(newden - logden) < TOL_1F1)) {
      logden = newden;
      return -exp(lognum - logden);
    }

    logden = newden;
    b++;

  }/*FOR*/

  /* large z: d/db log 1F1(1; b; z) -> digamma(b) - log(z). */
  return digamma(b0) - logz;

}/*DLOG1F1_DB*/
