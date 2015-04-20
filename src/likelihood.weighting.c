#include "include/rcore.h"
#include "include/loss.h"

void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    int debuglevel) {

int i = 0;
double maxw = 0;

  /* ensure the buffer is clean. */
  memset(w, '\0', n * sizeof(double));
  /* compute log-probabilities for each particle. */
  c_entropy_loss(fitted, data, n, TRUE, w, keep, FALSE, debuglevel);
  /* rescale before exponentiating them into probabilities. */
  maxw = w[d_which_max(w, n) - 1];

  for (i = 0; i < n; i++)
    w[i] = exp(w[i] - maxw);

}/*C_LW_WEIGHTS*/

SEXP lw_weights(SEXP fitted, SEXP data, SEXP keep, SEXP debug) {

int n = length(VECTOR_ELT(data, 0));
SEXP weights;

  PROTECT(weights = allocVector(REALSXP, n));

  c_lw_weights(fitted, data, n, REAL(weights), keep, isTRUE(debug));

  UNPROTECT(1);

  return weights;

}/*LW_WEIGHTS*/

