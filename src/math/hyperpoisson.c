#include "../include/rcore.h"
#include "../include/globals.h"
#include "./hypergeometric.h"
#include "./hyperpoisson.h"

/* log-rising factorial. */
#define lrf(z, n_val) (lgammafn(z + n_val) - lgammafn(z))

/* hard iteration cap for the inverse-CDF sampler, a last-resort backstop in
 * case the terms do not vanish as expected (e.g. non-finite parameters). */
#define MAXITER_RHYPOIS 10000000

/* density of a hyper-Poisson random variable. */
double dhypois(double x, double intensity, double dispersion, bool give_log) {

double logf = 0;

  if (ISNAN(x) || ISNAN(intensity) || ISNAN(dispersion))
    return NA_REAL;

  logf = -log1F1(1, dispersion, intensity);

  if (x > 0)
    logf += - lrf(dispersion, x) + x * log(intensity);

  if (give_log)
    return logf;
  else
    return exp(logf);

}/*DHYPOIS*/

/* mean and variance of a hyper-Poisson random variable. The recurrence
 * (dispersion + n) p_{n+1} = intensity p_n between successive probabilities
 * yields both moments in closed form from the normalising constant alone:
 * with p0 = 1 / 1F1(1; dispersion; intensity),
 *   mean     = intensity - (dispersion - 1) (1 - p0),
 *   variance = intensity - (dispersion - 1) p0 mean,
 * so a single 1F1 evaluation (which carries its own large-intensity asymptotic)
 * replaces the three moment series S_k = sum_n n^k intensity^n / (dispersion)_n.
 * Used by the EM/IRLS estimator to form the Fisher-scoring weights for the count
 * component. */
void hypois_moments(double intensity, double dispersion, double *mu,
    double *var) {

double p0 = 0, m = 0, v = 0;

  /* singular distribution with a spike at zero. */
  if (intensity < MACHINE_TOL) {

    *mu = 0;
    *var = 0;
    return;

  }/*THEN*/

  p0 = exp(-log1F1(1, dispersion, intensity));
  m = intensity - (dispersion - 1) * (1 - p0);
  v = intensity - (dispersion - 1) * p0 * m;

  *mu = m;
  *var = (v > 0) ? v : 0;

}/*HYPOIS_MOMENTS*/

/* mean of a hyper-Poisson random variable: lambda * d/dlambda log 1F1(1; b; l) =
 * (lambda / dispersion) * 1F1(2; dispersion + 1; lambda) / 1F1(1; dispersion;
 * lambda), evaluated in log-space. This is the cheaper mean-only form (two 1F1
 * evaluations) used for fitted values and predictions; the full moments are
 * available from hypois_moments(). */
double hypois_mean(double intensity, double dispersion) {

  return (intensity / dispersion) *
           exp(log1F1(2, dispersion + 1, intensity) -
               log1F1(1, dispersion, intensity));

}/*HYPOIS_MEAN*/

/* random sampling from a hyper-Poisson random variable. */
double rhypois(double intensity, double dispersion) {

int x = 0;
double log_u = 0, log_p = R_NegInf, log_pmf = 0, log_f11 = 0, log_eps = 0;

  if (ISNAN(intensity) || ISNAN(dispersion))
    return NA_REAL;

  /* singular distribution with a spike at zero. */
  if (intensity < MACHINE_TOL)
    return 0;

  /* dispersion ~= 1: standard Poisson. */
  if (fabs(dispersion - 1) < MACHINE_TOL)
    return rpois(intensity);

  /* dispersion ~= 0: the rising factorial (dispersion)_x becomes singular and
   * the hyper-Poisson tends to a shifted Poisson, X - 1 ~ Poisson(intensity);
   * use the closed form rather than the (ill-conditioned) inverse-CDF below. */
  if (dispersion < sqrt(MACHINE_TOL))
    return rpois(intensity) + 1;

  /* none of the shortcuts above need a uniform draw or the (expensive)
   * normalising constant, so only compute them here, for the general case. */
  log_u = log(unif_rand());
  log_f11 = log1F1(1, dispersion, intensity);
  log_eps = log(DBL_EPSILON);

  /* sample a uniform variate, walk through the log-CDF to identify the
   * corresponding quantile and the associated values. */
  while (log_p < log_u) {

    log_pmf = x * log(intensity) - log_f11 - lrf(dispersion, x);

    if (log_p == R_NegInf)
      log_p = log_pmf;
    else
      log_p = ((log_p > log_pmf) ? log_p : log_pmf) +
                log1p(exp(-fabs(log_p - log_pmf)));

    x++;

    /* the rising factorial makes the terms vanish super-exponentially, so the
     * series always converges; stop once the current term is too small to
     * change the accumulated probability (which also catches an imperfectly-
     * normalised tail when the sampled uniform is very close to one), or at the
     * hard iteration cap as a last resort. */
    if ((log_pmf - log_p < log_eps) || (x > MAXITER_RHYPOIS))
      break;

  }/*WHILE*/

  return x - 1;

}/*RHYPOIS*/
