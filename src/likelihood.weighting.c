#include "include/rcore.h"
#include "include/loss.h"

void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    int debuglevel) {

int i = 0, max_el = 0;
double maxw = 0;

  /* ensure the buffer is clean. */
  memset(w, '\0', n * sizeof(double));
  /* compute log-probabilities for each particle. */
  c_entropy_loss(fitted, data, n, TRUE, w, keep, FALSE, FALSE, debuglevel);
  /* rescale before exponentiating them into probabilities (if possible). */
  max_el = d_which_max(w, n);

  if (max_el == NA_INTEGER)
    return;
  else if ((max_el == 1) && (w[0] == R_NegInf))
    memset(w, '\0', n * sizeof(double));
  else {

    maxw = w[max_el - 1];
    for (i = 0; i < n; i++)
      w[i] = exp(w[i] - maxw);

  }/*ELSE*/

}/*C_LW_WEIGHTS*/

SEXP lw_weights(SEXP fitted, SEXP data, SEXP keep, SEXP debug) {

int n = length(VECTOR_ELT(data, 0));
SEXP weights;

  PROTECT(weights = allocVector(REALSXP, n));

  c_lw_weights(fitted, data, n, REAL(weights), keep, isTRUE(debug));

  UNPROTECT(1);

  return weights;

}/*LW_WEIGHTS*/

