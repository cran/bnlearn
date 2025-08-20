#include "../../include/rcore.h"
#include "../../include/globals.h"
#include "../../core/allocations.h"
#include "../../core/contingency.tables.h"
#include "../../math/linear.algebra.h"
#include "../parameters.h"

static long double estimate_loglik_tau(double *kappa, double tau, double alpha0,
    double s, int x_dim, int y_dim) {

long double loglik = 0;

  for (int i = 0; i < x_dim; i++)
    loglik += (digamma(tau * kappa[i]) - digamma(tau)) *
              (alpha0 - (tau * kappa[i]) + y_dim * (1 - s * kappa[i])) +
              lgammafn(tau * kappa[i]);

  loglik += - lgammafn(tau) - (y_dim * s * (x_dim - 1) / tau);

  return loglik;

}/*ESTIMATE_LOGLIK_TAU*/

static double estimate_loglik_kappa(double *kappa, double *digamma_rowsums,
    double tau, double s, double alpha0, int x_dim, int y_dim) {

long double loglik = 0;

  for (int i = 0; i < x_dim; i++) {

    loglik +=
      digamma(tau * kappa[i]) *
      (alpha0 - (tau * kappa[i]) + y_dim * (1 - s * kappa[i])) +
      lgammafn(tau * kappa[i]) -
      y_dim * (lgammafn(s * kappa[i]) + (1 - s * kappa[i]) * log(kappa[i])) +
      s * kappa[i] * digamma_rowsums[i];

  }/*FOR*/

  return (double) loglik;

}/*ESTIMATE_LOGLIK_KAPPA*/

static long double estimate_loglik_kappa_and_tau(double *nu, double *kappa,
    double tau, double alpha0, double s, int x_dim, int y_dim) {

double *digamma_rowsums = 0;
long double loglik_kappa_tau = 0;

  digamma_rowsums = (double *) Calloc1D(x_dim, sizeof(double));
  for (int i = 0; i < x_dim; i++)
    for (int j = 0; j < y_dim; j++)
      digamma_rowsums[i] += digamma(nu[CMC(i, j, x_dim)]);

  /* the joint log-likelihood of kappa and tau computed as the sum of the
   * individual log-likelihoods of kappa and tau, subtracting shared terms that
   * appear in both. */
  loglik_kappa_tau =
    estimate_loglik_kappa(kappa, digamma_rowsums, tau, s, alpha0, x_dim, y_dim) +
    estimate_loglik_tau(kappa, tau, alpha0, s, x_dim, y_dim);

  for (int i = 0; i < x_dim; i++)
    loglik_kappa_tau -=
      digamma(tau * kappa[i]) *
        (alpha0 - tau * kappa[i] - y_dim * (s * kappa[i] - 1)) +
       lgammafn(tau * kappa[i]);

  Free1D(digamma_rowsums);

  return loglik_kappa_tau;

}/*ESTIMATE_LOGLIK_KAPPA_AND_TAU*/

/* compute the log-likelihood component for nu. */
static long double estimate_loglik_nu(double *nu, double *kappa, cmcmap counts,
    double s) {

double colsum = 0;
long double loglik_nu = 0;

  for (int j = 0; j < counts.ncols; j++) {

    colsum = 0;
    for (int i = 0; i < counts.nrows; i++)
      colsum += nu[CMC(i, j, counts.nrows)];

    for (int i = 0; i < counts.nrows; i++)
      loglik_nu +=
        lgammafn(nu[CMC(i, j, counts.nrows)]) -
        (nu[CMC(i, j, counts.nrows)] - 1) * (digamma(nu[CMC(i, j, counts.nrows)]) - digamma(colsum)) -
        digamma(nu[CMC(i, j, counts.nrows)]) +
        CMEL(counts, i, j) * (digamma(nu[CMC(i, j, counts.nrows)]) - digamma(colsum)) +
        s * kappa[i] * digamma(nu[CMC(i, j, counts.nrows)]);

    loglik_nu += - lgammafn(colsum) - (s - counts.nrows) * digamma(colsum);

  }/*FOR*/

  return loglik_nu;

}/*ESTIMATE_LOGLIK_NU*/

static long double estimate_global_loglik(double *nu, double *kappa, double tau,
    double alpha0, double s, cmcmap counts) {

double *digamma_rowsums = NULL;
long double loglik = 0, loglik_kappa_and_tau = 0, loglik_nu = 0;

  digamma_rowsums = (double *) Calloc1D(counts.nrows, sizeof(double));
  for (int i = 0; i < counts.nrows; i++)
    for (int j = 0; j < counts.ncols; j++)
      digamma_rowsums[i] += digamma(nu[CMC(i, j, counts.nrows)]);

  loglik_kappa_and_tau =
    estimate_loglik_kappa_and_tau(nu, kappa, tau, alpha0, s, counts.nrows, counts.ncols);
  loglik_nu =
    estimate_loglik_nu(nu, kappa, counts, s);

  /* the global likelihood is the sum of its components... */
  loglik = loglik_kappa_and_tau + loglik_nu +
      /* ... plus terms that are not part of any component (disregarding
       * lgammafn(s0))... */
      counts.ncols * lgammafn(s) - counts.nrows * lgammafn(alpha0);

  /* ... minus shared terms that are duplicated across components. */
  for (int i = 0; i < counts.nrows; i++)
    loglik -= s * kappa[i] * digamma_rowsums[i];

  Free1D(digamma_rowsums);

return loglik;

}/*ESTIMATE_GLOBAL_LOGLIK*/

static void update_tau(double *kappa, int x_dim, int y_dim, double alpha0,
    double s, double *tau, hdstatus *err) {

int newton_iter = 0, max_iter = 200;
double gradient = 1, hessian = 0;
double tau_start = 100, tau_cur = 100;
double log_tau = log(tau_start);

  /* newton iterations: limit the maximum number of iterations to make sure not
   * to be stuck in an infinite loop regardless of how ill-conditioned the data
   * are. */
  for (newton_iter = 0; newton_iter < max_iter; newton_iter++) {

    gradient = hessian = 0;

    /* the partial derivative of the log-likelihood with respect to tau, the
     * first formula in Appendix B.3.1. */
    for (int i = 0; i < x_dim; i++)
      gradient +=
        (trigamma(tau_cur * kappa[i]) * kappa[i] - trigamma(tau_cur)) *
        (alpha0 - (tau_cur * kappa[i]) + y_dim * (1 - s * kappa[i]));

    gradient += y_dim * s * (x_dim - 1) / tau_cur / tau_cur;

    /* the partial second order derivative of the log-likelihood with respect to
     * tau, the second formula in Appendix B.3.1. */
    for (int i = 0; i < x_dim; i++)
      hessian +=
        (kappa[i] * kappa[i] * tetragamma(tau_cur * kappa[i]) - tetragamma(tau_cur)) *
        (alpha0 - (tau_cur * kappa[i]) - y_dim * (s * kappa[i] - 1)) -
        kappa[i] * kappa[i] * trigamma(tau_cur * kappa[i]);

    hessian += trigamma(tau_cur) -
                 2 * (s / tau_cur / tau_cur / tau_cur) * (y_dim * (x_dim - 1));

    /* updating tau according to the last two formulas in Appendix B.3.1. */
    log_tau += - gradient / (hessian * tau_cur + gradient);
    tau_cur = exp(log_tau);

    /* tau is the scale argument of the dirichlet distribution in the top layer
     * of the variational model in Equation 12, so it must be strictly positive:
     * restart newton's algorithm from further afar if it appears otherwise. */
    if (tau_cur < MACHINE_TOL || ISNAN(tau_cur)) {

      tau_start = tau_start * 10;
      tau_cur = tau_start;
      log_tau = log(tau_cur);

      (*err).tau_is_zero = TRUE;

    }/*THEN*/

    /* stop iterating if the (first order) partial derivative converged to zero
     * because that is where the maximum is. */
    if (fabs(gradient) < MACHINE_TOL)
      break;

  }/*FOR*/

  /* successful convergence if we did not reach the iteration limit. */
  if (newton_iter == max_iter)
    (*err).tau_convergence_fail = TRUE;

  /* save the final estimate of tau. */
  *tau = tau_cur;

}/*UPDATE_TAU*/

static void partial_derivatives_kappa(double *kappa, double *digamma_rowsums,
    double tau, double s, double alpha0, int x_dim, int y_dim, double *gradient,
    double *hessian) {

  for (int i = 0; i < x_dim; i++) {

    /* the partial derivative of the log-likelihood with respect to kappa,
     * which is the first equation in Appendix B.3.2. */
    gradient[i] =
      tau * trigamma(tau * kappa[i]) *
      (alpha0 - (tau * kappa[i]) - y_dim * (s * kappa[i] - 1)) -
      s * y_dim * (digamma(tau * kappa[i]) + digamma(s * kappa[i])
        - log(kappa[i]) - 1) -
      y_dim / kappa[i] + s * digamma_rowsums[i];

    /* the partial second derivative of the log-likelihood with respect to
     * kappa, which is the second equation in Appendix B.3.2. */
    hessian[i] =
      (tau * tau) * tetragamma(tau * kappa[i]) *
      (alpha0 - (tau * kappa[i]) - y_dim *(s * kappa[i] - 1)) -
      tau * trigamma(tau * kappa[i]) * (tau + 2 * s * y_dim) -
      s * s * y_dim * trigamma(s * kappa[i]) +
      s * y_dim / kappa[i] + y_dim / (kappa[i] * kappa[i]);

  }/*FOR*/

}/*PARTIAL_DERIVATIVES_KAPPA*/

static void update_kappa(double *nu, double tau, int x_dim, int y_dim,
    double alpha0, double s, double *kappa, hdstatus *err) {

int newton_iter = 0, max_iter = 200;
double loglik_kappa = 0, loglik_kappa_old = 0;
double *digamma_rowsums = NULL, *gradient = NULL, *hessian = NULL;
double *new_kappa = NULL, *delta_kappa = NULL;
double sqr_newton_decrement = 10 * MACHINE_TOL, step_size = 1;
double invhsum = 0, goverhsum =  0;
double expected_increase = 0, limit = 0;

  digamma_rowsums = (double *) Calloc1D(x_dim, sizeof(double));
  for (int i = 0; i < x_dim; i++)
    for (int j = 0; j < y_dim; j++)
      digamma_rowsums[i] += digamma(nu[CMC(i, j, x_dim)]);

  for (int i = 0; i < x_dim; i++)
    kappa[i] = 1.0 / x_dim;

  gradient = (double *) Calloc1D(x_dim, sizeof(double));
  hessian = (double *) Calloc1D(x_dim, sizeof(double));
  delta_kappa = (double *) Calloc1D(x_dim, sizeof(double));
  new_kappa = (double *) Calloc1D(x_dim, sizeof(double));

  loglik_kappa =
    estimate_loglik_kappa(kappa, digamma_rowsums, tau, s, alpha0, x_dim, y_dim);

  /* constrained newton method: g[] are the gradients, h[] are the diagonal
   * elements of the hessian from Appendix B.3.2. */
  for (newton_iter = 0; newton_iter < max_iter; newton_iter++) {

    /* reset the accumulators used to compute the newton step. */
    invhsum = goverhsum = 0;
    loglik_kappa_old = loglik_kappa;
    sqr_newton_decrement = expected_increase = 0;
    step_size = 1;

    partial_derivatives_kappa(kappa, digamma_rowsums, tau, s, alpha0, x_dim,
        y_dim, gradient, hessian);

    for (int i = 0; i < x_dim; i++) {

      invhsum += 1 / hessian[i];
      goverhsum += gradient[i] / hessian[i];

    }/*FOR*/

    /* the constrained newton step in the last formula of Appendix B.3.2. */
    for (int i = 0; i < x_dim; i++) {

      delta_kappa[i] = (goverhsum / invhsum - gradient[i]) / hessian[i];
      sqr_newton_decrement -= hessian[i] * delta_kappa[i] * delta_kappa[i];
      expected_increase += gradient[i] * delta_kappa[i];

      if (delta_kappa[i] < 0) {

        limit = (kappa[i] - 1e-10) / (-delta_kappa[i]);
        if (step_size > limit)
          step_size = limit;

      }/*THEN*/

    }/*FOR*/

    /* operate on a backup copy of kappa, so that we can always restore the
     * original estimates. */
    memcpy(new_kappa, kappa, x_dim * sizeof(double));

    while (loglik_kappa < loglik_kappa_old + 0.4 * step_size * expected_increase &&
            step_size > MACHINE_TOL) {

      for (int i = 0; i < x_dim; i++)
        new_kappa[i] = kappa[i] + step_size * delta_kappa[i];

      loglik_kappa = estimate_loglik_kappa(new_kappa, digamma_rowsums, tau, s,
                       alpha0, x_dim, y_dim);

      /* reduce the stepping by 10%. */
      step_size *= 0.9;

    }/*WHILE*/

    /* save the updated kappa values to make them available to the caller. */
    memcpy(kappa, new_kappa, x_dim * sizeof(double));

    /* the objective function did not change by an appreciable amount. */
    if (sqr_newton_decrement < 2 * MACHINE_TOL)
      break;
    /* the step size is so small that the estimates cannot change by any
       appreciable amount. */
    if (step_size < MACHINE_TOL)
      break;

  }/*FOR*/

  /* successful convergence if we did not reach the iteration limit. */
  if (newton_iter == max_iter)
    (*err).kappa_convergence_fail = TRUE;

  /* the kappa are Dirichlet paramaters, and therefore must be bounded away from
   * zero to meet regularity conditions. */
  for (int i = 0; i < x_dim; i++)
    kappa[i] = kappa[i] < MACHINE_TOL ? MACHINE_TOL : kappa[i];

  Free1D(digamma_rowsums);
  Free1D(gradient);
  Free1D(hessian);
  Free1D(delta_kappa);
  Free1D(new_kappa);

}/*UPDATE_KAPPA*/

static void update_kappa_and_tau(double *nu, double *kappa, int x_dim,
    int y_dim, double alpha0, double s, double *tau, hdstatus *err,
    bool debugging) {

int newton_iter = 0, max_iter = 200;
long double loglik_kappa_tau = 0, loglik_kappa_tau_old = 0, relative_difference = 1;

  /* tau is resetted to zero very time, instead of being iteratively updated. */
  *tau = 0;

  /* newton iterations: limit the maximum number of iterations to make sure not
   * to be stuck in an infinite loop regardless of how ill-conditioned the data
   * are. */
  for (newton_iter = 0; newton_iter < max_iter; newton_iter++) {

    if (debugging)
      Rprintf("    > updating kappa and tau, iteration %d: ", newton_iter + 1);

    /* reset the log-likelihood component for tau and kappa. */
    loglik_kappa_tau = 0;

    /* estimate tau and compute the corresponding log-likelihood component. */
    update_tau(kappa, x_dim, y_dim, alpha0, s, tau, err);
    estimate_loglik_tau(kappa, *tau, alpha0, s, x_dim, y_dim);

    /* estimate kappa and compute the corresponding log-likelihood component. */
    update_kappa(nu, *tau, x_dim, y_dim, alpha0, s, kappa, err);

    loglik_kappa_tau =
      estimate_loglik_kappa_and_tau(nu, kappa, *tau, alpha0, s, x_dim, y_dim);

    if (debugging)
      Rprintf("the log-likelihood is %lf.\n", (double)loglik_kappa_tau);

    relative_difference = (loglik_kappa_tau - loglik_kappa_tau_old) /
                            fabsl(loglik_kappa_tau_old);
    loglik_kappa_tau_old = loglik_kappa_tau;

    /* the log-likelihood is not increasing any longer, which suggests it
     * reached its peak. */
    if (fabsl(relative_difference) < MACHINE_TOL)
      break;

  }/*FOR*/

  /* successful convergence if we did not reach the iteration limit. */
  if (newton_iter == max_iter)
    (*err).kappa_tau_convergence_fail = TRUE;

}/*UPDATE_KAPPA_AND_TAU*/

static void update_nu(cmcmap counts, double *kappa, double s, double *nu,
    bool debugging) {

  if (debugging)
    Rprintf("    > updating nu.\n");

  /* compute the updated nu from the updated kappa, the imaginary sample size
   * and the data. */
  for (int i = 0; i < counts.nrows; i++)
    for (int j = 0; j < counts.ncols; j++)
      nu[CMC(i, j, counts.nrows)] = s * kappa[i] + CMEL(counts, i, j);

}/*UPDATE_NU*/

hdstatus c_hierarchical_dirichlet_parameters(cmcmap counts, double alpha0,
    double s, bool debugging, double *nu) {

int em_iter = 0, max_iter = 200;
double tau = 100;
double *kappa = NULL;
long double relative_difference = 0, loglik_old = 0, loglik_cur = 0;
hdstatus err = { 0 };

  /* allocate and initialise the parameters of the variational estimator. */
  for (int i = 0; i < counts.nrows * counts.ncols; i++)
    nu[i] = 1.0 / counts.nrows;

  kappa = (double *) Calloc1D(counts.nrows, sizeof(double));
  for (int i = 0; i < counts.nrows; i++)
    kappa[i] = 1.0 / counts.nrows;

  /* expectation-maximization iterations: limit the maximum number of iterations
   * to make sure not to be stuck in an infinite loop regardless of how
   * ill-conditioned the data are. */
  for (em_iter = 0; em_iter < max_iter; em_iter++) {

    if (debugging)
      Rprintf("  > iteration %d.\n", em_iter + 1);

    /* within each iteration, first update kappa and tau and compute the
     * associated log-likelihood component... */
    update_kappa_and_tau(nu, kappa, counts.nrows, counts.ncols, alpha0, s, &tau,
      &err, debugging);

    /* ... update nu and the associated log-likelihood component, using the
     * updated values of kappa... */
    update_nu(counts, kappa, s, nu, debugging);

    /* ... then add up the pieces of the log-likelihood. */
    loglik_cur = estimate_global_loglik(nu, kappa, tau, alpha0, s, counts);

    if (debugging)
      Rprintf("  > the log-likelihood is now %lf.\n", (double)loglik_cur);

    /* check for convergence using the relative difference between the
     * log-likelihood in the current and previous iteration. */
    loglik_old = (loglik_old == 0.0) ? 1e-10 : loglik_old;
    relative_difference = (loglik_cur - loglik_old) / fabsl(loglik_old);
    loglik_old = loglik_cur;

    /* the log-likelihood is not increasing any longer, which suggests it
     * reached its peak. */
    if (fabsl(relative_difference) < MACHINE_TOL)
      break;

  }/*FOR*/

  /* successful convergence if we did not reach the iteration limit. */
  if (em_iter == max_iter)
    err.outer_em_convergence_fail = TRUE;

  Free1D(kappa);

  return err;

}/*C_HIERARCHICAL_DIRICHLET_PARAMETERS*/

