#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../include/globals.h"
#include "./hypergeometric.h"
#include "./hyperpoisson.h"
#include "linear.algebra.h"
#include "reweighted.least.squares.h"

/* clamp for the linear predictor before exponentiating / inverse-logit. */
#define ETA_MAX 30.0
/* bounds and controls for the negative-binomial dispersion Newton update. */
#define THETA_MIN 1e-4
#define THETA_MAX 1e6
#define THETA_MAX_ITER 25
#define THETA_TOL 1e-8
/* bounds for the hyper-Poisson dispersion. */
#define DISP_MIN 1e-6
#define DISP_MAX 1e6

/* numerically stable inverse-logit, clamped strictly inside (0, 1). */
double irls_inv_logit(double eta) {

double p = 0;

  if (eta >= 0)
    p = 1 / (1 + exp(-eta));
  else
    p = exp(eta) / (1 + exp(eta));

  if (p < MACHINE_TOL)
    return MACHINE_TOL;
  else if (p > 1 - MACHINE_TOL)
    return 1 - MACHINE_TOL;
  else
    return p;

}/*IRLS_INV_LOGIT*/

/* linear predictor with the predictor clamped to [-ETA_MAX, ETA_MAX]. */
double irls_linpred(const double *X, int n, int p, const double *coef, int i) {

double eta = 0;

  /* NA coefficients mark aliased columns. */
  for (int k = 0; k < p; k++)
    if (!ISNAN(coef[k]))
      eta += X[i + k * n] * coef[k];

  if (eta > ETA_MAX)
    return ETA_MAX;
  else if (eta < -ETA_MAX)
    return -ETA_MAX;
  else
    return eta;

}/*IRLS_LINPRED*/

/* weighted least squares with just the intercept: the coefficient is the
 * weighted mean of the working response. */
static void c_wls0(const double *z, const double *w, int n, double *beta) {

long double sw = 0, swz = 0;

  for (int i = 0; i < n; i++) {

    sw += w[i];
    swz += w[i] * z[i];

  }/*FOR*/

  beta[0] = (double)(swz / sw);

}/*C_WLS0*/

/* weighted least squares with one regressor: closed-form weighted simple
 * regression. a constant regressor is collinear with the intercept, so it is
 * aliased to NA and the intercept absorbs the weighted mean. */
static void c_wls1(const double *X, const double *z, const double *w, int n,
    double *beta) {

const double *x1 = X + n;
long double sw = 0, swz = 0, sw1 = 0, s11 = 0, s1z = 0;
double zbar = 0, x1bar = 0, v11 = 0, c1 = 0;

  /* weighted means. */
  for (int i = 0; i < n; i++) {

    sw += w[i];
    swz += w[i] * z[i];
    sw1 += w[i] * x1[i];

  }/*FOR*/
  zbar = (double)(swz / sw);
  x1bar = (double)(sw1 / sw);

  /* weighted centred (co)variances. */
  for (int i = 0; i < n; i++) {

    double d1 = x1[i] - x1bar;
    s11 += w[i] * d1 * d1;
    s1z += w[i] * d1 * (z[i] - zbar);

  }/*FOR*/
  v11 = (double)(s11 / sw);
  c1 = (double)(s1z / sw);

  if (v11 < MACHINE_TOL) {

    beta[0] = zbar;
    beta[1] = NA_REAL;

  }/*THEN*/
  else {

    beta[1] = c1 / v11;
    beta[0] = zbar - beta[1] * x1bar;

  }/*ELSE*/

}/*C_WLS1*/

/* weighted least squares with two regressors: closed-form weighted estimates,
 * with the same three collinear configurations handled as in c_ols2(). */
static void c_wls2(const double *X, const double *z, const double *w, int n,
    double *beta) {

const double *x1 = X + n, *x2 = X + 2 * n;
long double sw = 0, swz = 0, sw1 = 0, sw2 = 0;
long double s11 = 0, s22 = 0, s12 = 0, s1z = 0, s2z = 0;
double zbar = 0, x1bar = 0, x2bar = 0;
double v11 = 0, v22 = 0, v12 = 0, c1 = 0, c2 = 0, den = 0;
int singular1 = FALSE, singular2 = FALSE;

  /* weighted means. */
  for (int i = 0; i < n; i++) {

    sw += w[i];
    swz += w[i] * z[i];
    sw1 += w[i] * x1[i];
    sw2 += w[i] * x2[i];

  }/*FOR*/
  zbar = (double)(swz / sw);
  x1bar = (double)(sw1 / sw);
  x2bar = (double)(sw2 / sw);

  /* weighted centred (co)variances. */
  for (int i = 0; i < n; i++) {

    double d1 = x1[i] - x1bar, d2 = x2[i] - x2bar, dz = z[i] - zbar;
    s11 += w[i] * d1 * d1;
    s22 += w[i] * d2 * d2;
    s12 += w[i] * d1 * d2;
    s1z += w[i] * d1 * dz;
    s2z += w[i] * d2 * dz;

  }/*FOR*/
  v11 = (double)(s11 / sw);
  v22 = (double)(s22 / sw);
  v12 = (double)(s12 / sw);
  c1 = (double)(s1z / sw);
  c2 = (double)(s2z / sw);

  /* there are three possible collinear configurations:
   *   1) the first regressor is constant and collinear with the intercept;
   *   2) the second regressor is constant, or collinear with the first;
   *   3) both. */
  singular1 = (v11 < MACHINE_TOL);
  singular2 = (v22 < MACHINE_TOL) ||
              (fabs(v12) / sqrt(v11 * v22) > 1 - MACHINE_TOL);

  if (singular1 && !singular2) {

    beta[1] = NA_REAL;
    beta[2] = c2 / v22;
    beta[0] = zbar - beta[2] * x2bar;

  }/*THEN*/
  else if (!singular1 && singular2) {

    beta[1] = c1 / v11;
    beta[2] = NA_REAL;
    beta[0] = zbar - beta[1] * x1bar;

  }/*THEN*/
  else if (singular1 && singular2) {

    beta[1] = NA_REAL;
    beta[2] = NA_REAL;
    beta[0] = zbar;

  }/*THEN*/
  else {

    den = v11 * v22 - v12 * v12;
    beta[1] = (v22 * c1 - v12 * c2) / den;
    beta[2] = (v11 * c2 - v12 * c1) / den;
    beta[0] = zbar - beta[1] * x1bar - beta[2] * x2bar;

  }/*ELSE*/

}/*C_WLS2*/

/* weighted least squares, special-casing small models with closed-form
 * estimates (as c_ols() does) and falling back on the pivoted QR otherwise. */
int c_weighted_least_squares(const double *X, int n, int p, const double *w,
    const double *z, double *beta, double *Xs, double *zs) {

double sd = 0;

  /* closed-form estimates for the null, one- and two-regressor models, avoiding
   * the QR decomposition (and its allocations) exactly as c_ols() does. */
  if (p == 1) {

    c_wls0(z, w, n, beta);
    return 0;

  }/*THEN*/
  else if (p == 2) {

    c_wls1(X, z, w, n, beta);
    return 0;

  }/*THEN*/
  else if (p == 3) {

    c_wls2(X, z, w, n, beta);
    return 0;

  }/*THEN*/

  /* row-scale Xs = sqrt(w) X and zs = sqrt(w) z, turning the weighted problem
   * into the equivalent OLS. */
  for (int i = 0; i < n; i++) {

    double sw = sqrt(w[i]);
    zs[i] = sw * z[i];
    for (int k = 0; k < p; k++)
      Xs[i + k * n] = sw * X[i + k * n];

  }/*FOR*/

  /* Xs is overwritten by its QR decomposition; fitted values and residuals are
   * not needed here. */
  c_qr(Xs, zs, n, p, NULL, NULL, beta, &sd);

  return 0;

}/*C_WEIGHTED_LEAST_SQUARES*/

/* weighted logistic regression by IRLS. */
int c_irls_logistic(const double *X, int n, int p, const double *y,
    const double *prior_w, double *beta, int max_iter, double tol,
    irls_scratch *ws) {

int iter = 0;
double *w = ws->w, *z = ws->z, *beta_new = ws->beta_new, *Xs = ws->Xs, *zs = ws->zs;

  /* an intercept-only logistic regression has the closed-form maximum likelihood
   * estimate beta0 = logit(weighted mean of the response), so skip the IRLS loop
   * entirely. */
  if (p == 1) {

    long double sw = 0, swy = 0;

    for (int i = 0; i < n; i++) {

      double pw = prior_w ? prior_w[i] : 1.0;
      sw += pw;
      swy += pw * y[i];

    }/*FOR*/

    double pbar = (double)(swy / sw);
    /* clamp strictly inside (0, 1) to guard against separation, as
     * irls_inv_logit() does inside the loop. */
    if (pbar < MACHINE_TOL)
      pbar = MACHINE_TOL;
    else if (pbar > 1 - MACHINE_TOL)
      pbar = 1 - MACHINE_TOL;
    beta[0] = log(pbar / (1 - pbar));

    return 1;

  }/*THEN*/

  for (iter = 0; iter < max_iter; iter++) {

    /* form the binomial working weights and responses at the current beta. */
    for (int i = 0; i < n; i++) {

      double eta = irls_linpred(X, n, p, beta, i);
      double pi = irls_inv_logit(eta);
      double vw = pi * (1 - pi);
      double pw = prior_w ? prior_w[i] : 1.0;

      w[i] = pw * vw;
      z[i] = eta + (y[i] - pi) / vw;

    }/*FOR*/

    c_weighted_least_squares(X, n, p, w, z, beta_new, Xs, zs);

    /* an aliased column stays NA in both beta_new and beta, so the NaN gap it
     * produces is ignored by the comparison below (NaN > delta is false). */
    double delta = 0;
    for (int k = 0; k < p; k++) {

      double d = fabs(beta_new[k] - beta[k]);
      if (d > delta)
        delta = d;
      beta[k] = beta_new[k];

    }/*FOR*/

    if (delta < tol) {

      iter++;
      break;

    }/*THEN*/

  }/*FOR*/

  return iter;

}/*C_IRLS_LOGISTIC*/

/* Newton update of the negative-binomial dispersion theta. */
static double nb_update_theta(const double *X, int n, int p, const double *y,
    const double *prior_w, const double *beta, double th) {

  for (int it = 0; it < THETA_MAX_ITER; it++) {

    double score = 0, info = 0;

    for (int i = 0; i < n; i++) {

      double pw = prior_w ? prior_w[i] : 1.0;
      if (pw <= 0)
        continue;

      double mu = exp(irls_linpred(X, n, p, beta, i));
      if (mu < MACHINE_TOL)
        mu = MACHINE_TOL;
      double tpm = th + mu;

      score += pw * (digamma(th + y[i]) - digamma(th) + log(th) + 1
                     - log(tpm) - (y[i] + th) / tpm);
      info  += pw * (- trigamma(th + y[i]) + trigamma(th) - 1.0 / th
                     + 2.0 / tpm - (y[i] + th) / (tpm * tpm));

    }/*FOR*/

    /* a non-positive or non-finite information makes the Newton step
     * unreliable, so stop and keep the current theta. */
    if (!R_FINITE(score) || !R_FINITE(info) || info < MACHINE_TOL)
      break;

    double th_new = th + score / info;
    if (th_new < THETA_MIN)
      th_new = THETA_MIN;
    else if (th_new > THETA_MAX)
      th_new = THETA_MAX;

    if (fabs(th_new - th) < THETA_TOL * th) {

      th = th_new;
      break;

    }/*THEN*/

    th = th_new;

  }/*FOR*/

  return th;

}/*NB_UPDATE_THETA*/

/* weighted profile hyper-Poisson score d/d(dispersion). */
static double hp_disp_score(const double *X, int n, int p, const double *y,
    const double *prior_w, const double *beta, double disp) {

double score = 0;

  for (int i = 0; i < n; i++) {

    double w = prior_w ? prior_w[i] : 1.0;
    if (w <= 0)
      continue;
    double lambda = exp(irls_linpred(X, n, p, beta, i));
    if (lambda < MACHINE_TOL)
      lambda = MACHINE_TOL;
    score += w * (-(digamma(disp + y[i]) - digamma(disp)) -
                  dlog1F1_db(disp, lambda));

  }/*FOR*/

  return score;

}/*HP_DISP_SCORE*/

/* Newton step for the hyper-Poisson dispersion on the weighted
 * profile score; the dispersion is kept within [DISP_MIN, DISP_MAX]. */
static double hp_update_dispersion(const double *X, int n, int p,
    const double *y, const double *prior_w, const double *beta, double disp) {

double s0 = 0, db = 0, s1 = 0, deriv = 0;

  s0 = hp_disp_score(X, n, p, y, prior_w, beta, disp);
  db = disp * 1e-4 + 1e-8;
  s1 = hp_disp_score(X, n, p, y, prior_w, beta, disp + db);
  deriv = (s1 - s0) / db;

  if (!R_FINITE(s0) || !R_FINITE(deriv) || fabs(deriv) < MACHINE_TOL)
    return disp;

  double disp_new = disp - s0 / deriv;
  if (disp_new < DISP_MIN)
    return DISP_MIN;
  else if (disp_new > DISP_MAX)
    return DISP_MAX;
  else
    return disp_new;

}/*HP_UPDATE_DISPERSION*/

/* weighted count regression by Fisher-scoring IRLS, jointly updating the
 * dispersion, for the family selected by `family`. */
int c_irls_count(glm_family_e family, const double *X, int n, int p,
    const double *y, const double *prior_w, double *beta, double *dispersion,
    int max_iter, double tol, irls_scratch *ws) {

int iter = 0;
double disp = *dispersion;
double *w = ws->w, *z = ws->z, *beta_new = ws->beta_new, *Xs = ws->Xs, *zs = ws->zs;

  for (iter = 0; iter < max_iter; iter++) {

    /* form the family-specific working weights and responses (log link). */
    for (int i = 0; i < n; i++) {

      double eta = irls_linpred(X, n, p, beta, i);
      double pw = prior_w ? prior_w[i] : 1.0;
      double mu = 0, weight = 0, deriv = 0;

      if (family == GLM_NEGBIN) {

        /* mean exp(eta); dmu/deta = mu; weight = mu * theta / (theta + mu). */
        mu = exp(eta);
        if (mu < MACHINE_TOL)
          mu = MACHINE_TOL;
        weight = mu * disp / (disp + mu);
        deriv = mu;

      }/*THEN*/
      else if (family == GLM_HYPERPOISSON) {

        /* dmu/deta = variance (one-parameter exponential family in
         * eta = log(lambda)), so the weight equals the variance. */
        double lambda = exp(eta), var = 0;
        if (lambda < MACHINE_TOL)
          lambda = MACHINE_TOL;
        hypois_moments(lambda, disp, &mu, &var);
        if (var < MACHINE_TOL)
          var = MACHINE_TOL;
        weight = var;
        deriv = var;

      }/*THEN*/

      w[i] = pw * weight;
      z[i] = eta + (y[i] - mu) / deriv;

    }/*FOR*/

    c_weighted_least_squares(X, n, p, w, z, beta_new, Xs, zs);

    double delta = 0;
    for (int k = 0; k < p; k++) {

      double d = fabs(beta_new[k] - beta[k]);
      if (d > delta)
        delta = d;
      beta[k] = beta_new[k];

    }/*FOR*/

    /* update the dispersion at the new coefficients. */
    if (family == GLM_NEGBIN)
      disp = nb_update_theta(X, n, p, y, prior_w, beta, disp);
    else if (family == GLM_HYPERPOISSON)
      disp = hp_update_dispersion(X, n, p, y, prior_w, beta, disp);

    if (delta < tol) {

      iter++;
      break;

    }/*THEN*/

  }/*FOR*/

  *dispersion = disp;

  return iter;

}/*C_IRLS_COUNT*/
