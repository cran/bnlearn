#include "../include/rcore.h"
#include "../include/globals.h"
#include "../core/math.functions.h"
#include "../fitted/fitted.h"
#include "../core/data.table.h"
#include "loglikelihood/loglikelihood.h"
#include "../minimal/common.h"

void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    bool debugging) {

int i = 0, max_el = 0;
double maxw = 0;
fitted_bn bn = fitted_network_from_SEXP(fitted);
SEXP metadata, complete_nodes, nodes_in_fitted, keep2;

  memset(w, '\0', n * sizeof(double));
  /* find out which nodes to compute the log-likelihood for, zero-indexed. */
  PROTECT(nodes_in_fitted = getAttrib(fitted, R_NamesSymbol));
  PROTECT(keep2 = match(keep, nodes_in_fitted, 0));

  /* extract the metadata from the data. */
  PROTECT(metadata = getAttrib(data, BN_MetaDataSymbol));
  PROTECT(complete_nodes = getListElement(metadata, "complete.nodes"));

  if (bn.type == DNET || bn.type == ONET || bn.type == DONET) {

    ddata dt = ddata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep2);

    bysample_discrete_loglikelihood(bn, dt, w, debugging);

  }/*THEN*/
  else if (bn.type == GNET) {

    gdata dt = gdata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep2);

    bysample_gaussian_loglikelihood(bn, dt, w, TRUE, debugging);

  }/*THEN*/
  else if (bn.type == CGNET) {

    cgdata dt = cgdata_from_SEXP(data, 0, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep2);

    bysample_clgaussian_loglikelihood(bn, dt, w, TRUE, debugging);

  }/*THEN*/

  FreeFittedBN(bn);

  UNPROTECT(4);

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

