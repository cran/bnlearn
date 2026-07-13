#include "../../include/rcore.h"
#include "../../core/data.table.h"
#include "../../math/sparse.amat.h"
#include "../../minimal/common.h"
#include "../../parameters/parameters.h"
#include "../graphs.h"

/* get the number of parameters of the whole network (mixed case, also handles
 * discrete and Gaussian networks). */
SEXP nparams_structure(SEXP graph, SEXP data, SEXP estimator, SEXP debug) {

double node_params = 0, all_params = 0;
bool debugging = isTRUE(debug);
estimator_e est = estimator_to_enum(CHAR(STRING_ELT(estimator, 0)));
tabular dt = tabular_from_SEXP(data, 0, 0);
/* the transposed adjacency matrix gives each node's parents (as node indexes)
 * directly, without matching the parent labels one node at a time. */
sparse_amat parents = sparse_amat_from_SEXP(graph, TRUE);

  for (int i = 0; i < parents.dim; i++) {

    node_params = nparams_node_count(parents, dt, i, est);

    if (debugging)
      Rprintf("* node %s has %.0lf parameter(s).\n", parents.labels[i],
        node_params);

    all_params += node_params;

  }/*FOR*/

  FreeTAB(dt);
  FreeSparseAMAT(parents);

  return ScalarReal(all_params);

}/*NPARAMS_STRUCTURE*/

