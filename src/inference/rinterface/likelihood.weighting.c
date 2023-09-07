#include "../../include/rcore.h"
#include "../../include/sampling.h"

/* compute the weights of the likelihood weighting particles. */
SEXP lw_weights(SEXP fitted, SEXP data, SEXP keep, SEXP debug) {

int n = length(VECTOR_ELT(data, 0));
SEXP weights;

  PROTECT(weights = allocVector(REALSXP, n));

  c_lw_weights(fitted, data, n, REAL(weights), keep, isTRUE(debug));

  UNPROTECT(1);

  return weights;

}/*LW_WEIGHTS*/

