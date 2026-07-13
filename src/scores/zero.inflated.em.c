#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/data.table.h"
#include "../include/globals.h"
#include "../math/hyperpoisson.h"
#include "../math/reweighted.least.squares.h"
#include "../minimal/common.h"
#include "../minimal/data.frame.h"
#include "scores.h"

/* inner IRLS controls for the full M-step (the EM controls come from R). */
#define IRLS_MAX_ITER 50
#define IRLS_TOL 1e-8
/* maximum step-halvings for the one-step (generalized EM) M-step safeguard. */
#define MAX_HALVING 10
/* floor in the relative parameter-change stopping rule. */
#define PAR_FLOOR 1e-3

/* log of the count density for the family: negative binomial (mean m = mu,
 * dispersion = theta) or hyper-Poisson (mean argument m = lambda, dispersion).
 * used by the E-step zero probability and the observed-data log-likelihood. */
static inline double count_log_density(glm_family_e family, double y, double m,
    double disp) {

  if (family == GLM_NEGBIN)
    return lgammafn(y + disp) - lgammafn(disp) - lgammafn(y + 1) +
             disp * (log(disp) - log(disp + m)) + y * (log(m) - log(disp + m));
  else
    return dhypois(y, m, disp, TRUE);

}/*COUNT_LOG_DENSITY*/

/* observed-data log-likelihood of the zero-inflated count distribution. */
static double em_count_loglik(glm_family_e family, const double *X, int n, int p,
    const double *y, const double *gamma, const double *beta, double disp) {

long double llik = 0;

  for (int i = 0; i < n; i++) {

    double pi = irls_inv_logit(irls_linpred(X, n, p, gamma, i));
    double m = exp(irls_linpred(X, n, p, beta, i));
    if (m < MACHINE_TOL)
      m = MACHINE_TOL;

    if (y[i] == 0)
      llik += log(pi + (1 - pi) * exp(count_log_density(family, 0, m, disp)));
    else
      llik += log(1 - pi) + count_log_density(family, y[i], m, disp);

  }/*FOR*/

  return (double)llik;

}/*EM_COUNT_LOGLIK*/

/* starting values: an intercept-only working fit for the count mean, dispersion
 * = 1, and the excess of observed zeroes over the count-distribution expectation
 * for the zero-inflation intercept. all slopes start at zero. */
static void em_init(glm_family_e family, int n, int p, const double *y,
    double *gamma, double *beta, double *disp) {

double ybar = 0;
int nzero = 0;

  for (int i = 0; i < n; i++) {

    ybar += y[i];
    if (y[i] == 0)
      nzero++;

  }/*FOR*/
  ybar /= n;

  for (int k = 0; k < p; k++) {

    beta[k] = 0;
    gamma[k] = 0;

  }/*FOR*/

  double d = 1.0;
  double m0 = (ybar > 0.1) ? ybar : 0.1;
  beta[0] = log(m0);

  double p0_obs = (double)nzero / n;
  double p0_count = exp(count_log_density(family, 0, m0, d));
  double excess = p0_obs - p0_count;
  if (excess < 0.02)
    excess = 0.02;
  else if (excess > 0.95)
    excess = 0.95;
  gamma[0] = log(excess / (1 - excess));

  *disp = d;

}/*EM_INIT*/

/* write the internal (gamma, beta, dispersion) into the bnlearn estimates[]
 * layout. for the negative binomial the count model is parameterised by the
 * logit of the success probability, logit(p) = X.beta - log(theta), so the
 * intercept absorbs -log(theta); for the hyper-Poisson the intensity coefficients
 * are beta directly. both store the dispersion on the log scale. */
static void to_estimates(glm_family_e family, const double *gamma,
    const double *beta, double disp, double *estimates, int p) {

  for (int k = 0; k < p; k++)
    estimates[k] = gamma[k];

  if (family == GLM_NEGBIN) {

    estimates[p] = beta[0] - log(disp);
    for (int k = 1; k < p; k++)
      estimates[p + k] = beta[k];

  }/*THEN*/
  else {

    for (int k = 0; k < p; k++)
      estimates[p + k] = beta[k];

  }/*ELSE*/

  estimates[2 * p] = log(disp);

}/*TO_ESTIMATES*/

/* parameter estimation and local log-likelihood for a zero-inflated count node
 * by EM with the GLM components fitted by IRLS, shared by the negative-binomial
 * and hyper-Poisson families. it is the sole estimator for these node types,
 * both for parameter learning (mle-zihp, mle-zinb) and for structure learning
 * (the loglik/aic/bic/nal/pnal count scores); it fills the estimates[] layout
 * [ inflation (np+1) | count coefs (np+1) | log-dispersion ], np = 2 * nparents
 * + 3, that the rest of the package expects. with one_step = TRUE the M-step
 * takes a single (step-halving-safeguarded) IRLS step per component instead of
 * re-fitting each GLM to convergence. */
double em_irls_node(SEXP target, SEXP x, SEXP data, double *estimates,
    double *nparams, bool pnal, double k, bool debugging, int em_max_iter,
    double em_tol, bool one_step, glm_family_e family) {

double loglik = 0;
double *y = NULL, *ymiss = NULL, *fit_y = NULL;
double *X = NULL, *gamma = NULL, *beta = NULL, *tau = NULL, *pw = NULL;
double *old_gamma = NULL, *old_beta = NULL;
int np = 0, p = 0, n = 0;
bool use_sub = FALSE;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, local_data;
tabular dt = { 0 }, sub = { 0 }, fit_dt = { 0 };
irls_scratch ws = { 0 };

  /* get the node and its parents from the network. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  parents = getListElement(node_t, "parents");

  np = 2 * length(parents) + 3;   /* full parameter count. */
  p = length(parents) + 1;        /* coefficients per GLM component. */

  /* extract the response and the parents. */
  y = REAL(c_dataframe_column(data, target, TRUE, FALSE));
  PROTECT(local_data = c_dataframe_column(data, parents, FALSE, TRUE));
  dt = tabular_from_SEXP(local_data, 0, 0);
  fit_dt = dt;
  fit_y = y;

  /* subset locally-complete observations when there are missing values. */
  if (pnal) {

    sub = new_tabular(dt.m.nobs, 0, dt.m.ncols);
    bool *missing = Calloc1D(dt.m.nobs, sizeof(bool));
    tabular_incomplete_cases(&dt, missing, 0, 0);
    for (int i = 0; i < dt.m.nobs; i++)
      missing[i] |= ISNAN(y[i]);
    ymiss = Calloc1D(dt.m.nobs, sizeof(double));
    for (int i = 0, j = 0; i < dt.m.nobs; i++)
      if (!missing[i])
        ymiss[j++] = y[i];
    tabular_subsample_by_logical(&dt, &sub, missing, 0, 0);
    Free1D(missing);
    use_sub = TRUE;
    fit_dt = sub;
    fit_y = ymiss;

  }/*THEN*/

  n = fit_dt.m.nobs;

  if (n == 0) {

    /* unidentifiable: no locally-complete data. */
    if (estimates)
      for (int i = 0; i < np; i++)
        estimates[i] = NA_REAL;
    loglik = pnal ? NA_REAL : R_NegInf;

  }/*THEN*/
  else {

    /* design matrix X (column-major): intercept column then parents. */
    X = Calloc1D(n * p, sizeof(double));
    for (int i = 0; i < n; i++)
      X[i] = 1.0;
    for (int j = 0; j < fit_dt.m.ncols; j++)
      for (int i = 0; i < n; i++)
        X[i + (j + 1) * n] = fit_dt.ccol[j][i];

    gamma = Calloc1D(p, sizeof(double));
    beta = Calloc1D(p, sizeof(double));
    old_gamma = Calloc1D(p, sizeof(double));
    old_beta = Calloc1D(p, sizeof(double));
    tau = Calloc1D(n, sizeof(double));
    pw = Calloc1D(n, sizeof(double));

    /* IRLS scratch, allocated once and shared by both M-step components across
     * all EM iterations instead of being reallocated inside each IRLS call. */
    ws.w = Calloc1D(n, sizeof(double));
    ws.z = Calloc1D(n, sizeof(double));
    ws.beta_new = Calloc1D(p, sizeof(double));
    ws.Xs = Calloc1D(n * p, sizeof(double));
    ws.zs = Calloc1D(n, sizeof(double));

    double disp = 1.0;

    em_init(family, n, p, fit_y, gamma, beta, &disp);
    loglik = em_count_loglik(family, X, n, p, fit_y, gamma, beta, disp);

    /* EM outer loop. */
    for (int iter = 0; iter < em_max_iter; iter++) {

      /* snapshot the current parameters to monitor their relative change. */
      double old_disp = disp;
      memcpy(old_gamma, gamma, p * sizeof(double));
      memcpy(old_beta, beta, p * sizeof(double));

      /* E-step: posterior probability that a zero is a structural zero. */
      for (int i = 0; i < n; i++) {

        if (fit_y[i] == 0) {

          double pi = irls_inv_logit(irls_linpred(X, n, p, gamma, i));
          double m = exp(irls_linpred(X, n, p, beta, i));
          if (m < MACHINE_TOL)
            m = MACHINE_TOL;
          double p0 = exp(count_log_density(family, 0, m, disp));
          tau[i] = pi / (pi + (1 - pi) * p0);
          pw[i] = 1 - tau[i];

        }/*THEN*/
        else {

          tau[i] = 0;
          pw[i] = 1;

        }/*ELSE*/

      }/*FOR*/

      /* M-step: weighted logistic GLM (zero) and weighted count GLM. */
      double ll_new = 0;

      if (one_step) {

        /* generalized EM: a single IRLS step per component, with step-halving
         * on the observed-data log-likelihood to keep the EM monotone. */
        c_irls_logistic(X, n, p, tau, NULL, gamma, 1, IRLS_TOL, &ws);
        c_irls_count(family, X, n, p, fit_y, pw, beta, &disp, 1, IRLS_TOL, &ws);

        ll_new = em_count_loglik(family, X, n, p, fit_y, gamma, beta, disp);

        for (int h = 0; (ll_new < loglik) && (h < MAX_HALVING); h++) {

          for (int kk = 0; kk < p; kk++) {

            gamma[kk] = 0.5 * (gamma[kk] + old_gamma[kk]);
            beta[kk] = 0.5 * (beta[kk] + old_beta[kk]);

          }/*FOR*/
          disp = sqrt(disp * old_disp);
          ll_new = em_count_loglik(family, X, n, p, fit_y, gamma, beta, disp);

        }/*FOR*/

      }/*THEN*/
      else {

        /* full M-step: re-fit each weighted GLM to convergence. */
        c_irls_logistic(X, n, p, tau, NULL, gamma, IRLS_MAX_ITER, IRLS_TOL, &ws);
        c_irls_count(family, X, n, p, fit_y, pw, beta, &disp, IRLS_MAX_ITER, IRLS_TOL, &ws);

        ll_new = em_count_loglik(family, X, n, p, fit_y, gamma, beta, disp);

      }/*ELSE*/

      /* relative observed-data log-likelihood change... */
      double crit = fabs(ll_new - loglik) / (fabs(loglik) + 1.0);
      loglik = ll_new;

      /* ... and the largest relative parameter change (scale-invariant, floored
       * so near-zero coefficients do not dominate and so a coefficient drifting
       * to a large magnitude on a non-identified direction can still settle). */
      double dpar = fabs(disp - old_disp) / (fabs(old_disp) + PAR_FLOOR);
      for (int kk = 0; kk < p; kk++) {

        double dg = fabs(gamma[kk] - old_gamma[kk]) / (fabs(old_gamma[kk]) + PAR_FLOOR);
        double db = fabs(beta[kk] - old_beta[kk]) / (fabs(old_beta[kk]) + PAR_FLOOR);
        if (dg > dpar)
          dpar = dg;
        if (db > dpar)
          dpar = db;

      }/*FOR*/

      if (debugging)
        Rprintf("  > EM iteration %d, log-likelihood %lf.\n", iter + 1, loglik);

      if (crit < em_tol && dpar < em_tol)
        break;

    }/*FOR*/

    /* store the parameters in the bnlearn estimates[] layout. */
    if (estimates)
      to_estimates(family, gamma, beta, disp, estimates, p);

    /* node-average penalised likelihood for (future) structure-learning use. */
    if (pnal) {

      loglik /= n;
      loglik -= k / dt.m.nobs * np;

    }/*THEN*/

  }/*ELSE*/

  if (debugging)
    Rprintf("  > log-likelihood is %lf.\n", loglik);

  if (nparams)
    *nparams = np;

  /* cleanup. */
  if (X)
    Free1D(X);
  if (gamma)
    Free1D(gamma);
  if (beta)
    Free1D(beta);
  if (old_gamma)
    Free1D(old_gamma);
  if (old_beta)
    Free1D(old_beta);
  if (tau)
    Free1D(tau);
  if (pw)
    Free1D(pw);
  if (ws.w)
    Free1D(ws.w);
  if (ws.z)
    Free1D(ws.z);
  if (ws.beta_new)
    Free1D(ws.beta_new);
  if (ws.Xs)
    Free1D(ws.Xs);
  if (ws.zs)
    Free1D(ws.zs);
  FreeTAB(dt);
  if (use_sub) {

    FreeTAB(sub);
    Free1D(ymiss);

  }/*THEN*/

  UNPROTECT(1);

  return loglik;

}/*EM_IRLS_NODE*/
