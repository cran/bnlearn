#include "../../include/rcore.h"
#include "../../fitted/fitted.h"

/* compute the number of parameters of a fitted model. */
SEXP nparams_fitted(SEXP fitted, SEXP effective, SEXP debug) {

double node_params = 0, all_params = 0;
fitted_bn bn = fitted_network_from_SEXP(fitted);

  for (int i = 0; i < bn.nnodes; i++) {

    node_params = nparams_fitted_node(bn.ldists[i], bn.node_types[i],
                    isTRUE(effective));

    if (isTRUE(debug))
      Rprintf("* node %s has %.0lf parameter(s).\n", bn.labels[i], node_params);

    all_params += node_params;

  }/*FOR*/

  FreeFittedBN(bn);

  return ScalarReal(all_params);

}/*NPARAMS_FITTED*/

