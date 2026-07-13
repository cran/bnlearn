#include "../include/rcore.h"
#include "../core/data.table.h"
#include "../include/globals.h"
#include "scores.h"

/* translate the regression coefficients of a zero-inflated count node into the
 * canonical parameters of its local distribution. these transforms are shared by
 * parameter estimation, prediction, random generation and the log-likelihood, so
 * they live apart from any single estimator. the zero-inflation probability and
 * (for the negative binomial) the success probability are logistic transforms of
 * the linear predictors, clamped strictly inside (0, 1) so that downstream
 * likelihoods stay finite under separation; the hyper-Poisson intensity is the
 * exponential of its linear predictor. */

/* numerically stable logistic transform of a linear predictor, in place, clamped
 * strictly inside (0, 1). */
static void logistic_in_place(double *x, int nobs) {

  for (int i = 0; i < nobs; i++) {

    double prob = (x[i] >= 0) ? 1 / (1 + exp(-x[i])) :
                    exp(x[i]) / (1 + exp(x[i]));
    if (prob < MACHINE_TOL)
      x[i] = MACHINE_TOL;
    else if (prob > 1 - MACHINE_TOL)
      x[i] = 1 - MACHINE_TOL;
    else
      x[i] = prob;

  }/*FOR*/

}/*LOGISTIC_IN_PLACE*/

/* accumulate a linear predictor intercept + parents . coefficients, in place.
 * an NA coefficient marks a rank-deficient (aliased) parent and contributes
 * nothing, matching the way glm() drops such terms. */
static void linear_predictor(double *x, double *coefs, tabular data) {

int nobs = data.m.nobs;

  for (int i = 0; i < nobs; i++)
    x[i] = ISNAN(coefs[0]) ? 0 : coefs[0];
  for (int j = 0; j < data.m.ncols; j++)
    if (!ISNAN(coefs[j + 1]))
      for (int i = 0; i < nobs; i++)
        x[i] += data.ccol[j][i] * coefs[j + 1];

}/*LINEAR_PREDICTOR*/

/* zero-inflated hyper-Poisson: logistic zero-inflation probability, exponential
 * intensity. */
void zihp_coefs_to_params(double *zinf_prob, double *zinf_coefs,
    double *intensity, double *inten_coefs, tabular data) {

int nobs = data.m.nobs;

  linear_predictor(zinf_prob, zinf_coefs, data);
  logistic_in_place(zinf_prob, nobs);

  linear_predictor(intensity, inten_coefs, data);
  for (int i = 0; i < nobs; i++)
    intensity[i] = exp(intensity[i]);

}/*ZIHP_COEFS_TO_PARAMS*/

/* zero-inflated negative binomial: logistic zero-inflation and success
 * probabilities. */
void zinb_coefs_to_params(double *zinf_prob, double *zinf_coefs,
    double *succ_prob, double *succ_coefs, tabular data) {

int nobs = data.m.nobs;

  linear_predictor(zinf_prob, zinf_coefs, data);
  logistic_in_place(zinf_prob, nobs);

  linear_predictor(succ_prob, succ_coefs, data);
  logistic_in_place(succ_prob, nobs);

}/*ZINB_COEFS_TO_PARAMS*/
