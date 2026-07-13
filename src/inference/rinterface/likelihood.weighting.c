#include "../../include/rcore.h"
#include "../../core/data.table.h"
#include "../../fitted/fitted.h"
#include "../../include/globals.h"
#include "../../include/sampling.h"
#include "../../minimal/common.h"

/* build the tabular passed to c_lw_weights() from the SEXP data. */
tabular lw_data_from_SEXP(SEXP data, SEXP fitted, SEXP keep) {

tabular dt = tabular_from_SEXP(data, 0, 0);
SEXP keep2, complete;

  PROTECT(keep2 = match(keep, getAttrib(fitted, R_NamesSymbol), 0));
  complete = getListElement(getAttrib(data, BN_MetaDataSymbol), "complete.nodes");
  meta_copy_names(&(dt.m), 0, data);
  meta_init_flags(&(dt.m), 0, complete, keep2);
  UNPROTECT(1);

  return dt;

}/*LW_DATA_FROM_SEXP*/

/* compute the weights of the likelihood weighting particles. */
SEXP lw_weights(SEXP fitted, SEXP data, SEXP keep, SEXP debug) {

int n = length(VECTOR_ELT(data, 0));
fitted_bn bn = fitted_network_from_SEXP(fitted);
tabular dt = lw_data_from_SEXP(data, fitted, keep);
SEXP weights;

  PROTECT(weights = allocVector(REALSXP, n));

  c_lw_weights(bn, dt, REAL(weights), isTRUE(debug));

  FreeTAB(dt);
  FreeFittedBN(bn);

  UNPROTECT(1);

  return weights;

}/*LW_WEIGHTS*/
