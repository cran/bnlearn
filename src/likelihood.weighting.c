#include "common.h"

SEXP lw_weights(SEXP fitted, SEXP orig_data, SEXP by_sample, SEXP keep,
    SEXP debug) {

int i = 0, n = 0;
double maxw = 0, *w = NULL;
SEXP weights;

  /* compute log-probabilities for each particle. */
  PROTECT(weights = entropy_loss(fitted, orig_data, by_sample, keep, debug));
  n = length(weights);
  w = REAL(weights);
  /* resscale before exponentiating them into probabilities. */
  maxw = w[which_max(w, n) - 1];

  for (i = 0; i < n; i++)
    w[i] = exp(w[i] - maxw);

  UNPROTECT(1);

  return weights;

}/*LW_WEIGHTS*/

