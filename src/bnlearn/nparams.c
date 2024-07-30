#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/common.h"
#include "../math/linear.algebra.h"

/* get the number of parameters of the whole network (mixed case, also handles
 * discrete and Gaussian networks). */
SEXP nparams_cgnet(SEXP graph, SEXP data, SEXP debug) {

int i = 0, j = 0, nnodes = 0;
int *nlevels = NULL, *index = NULL, ngp = 0;
double nconfig = 0, node_params = 0, all_params = 0;
bool debugging = isTRUE(debug);
SEXP nodes = R_NilValue, node_data, parents, temp;

  /* get nodes' number and data. */
  node_data = getListElement(graph, "nodes");
  nnodes = length(node_data);
  PROTECT(nodes = getAttrib(node_data, R_NamesSymbol));
  /* cache the number of levels of each variables (zero = continuous). */
  nlevels = Calloc1D(nnodes, sizeof(int));
  for (i = 0; i < nnodes; i++) {

    temp = VECTOR_ELT(data, i);
    if (TYPEOF(temp) == INTSXP)
      nlevels[i] = NLEVELS(temp);

  }/*FOR*/

  for (i = 0; i < nnodes; i++) {

    /* extract the parents of the node and match them. */
    parents = getListElement(VECTOR_ELT(node_data, i), "parents");
    PROTECT(temp = match(nodes, parents, 0));
    index = INTEGER(temp);

    /* compute the number of regressors and of configurations. */
    for (j = 0, ngp = 0, nconfig = 1; j < length(parents); j++) {

      if (nlevels[index[j] - 1] == 0)
        ngp++;
      else
        nconfig *= nlevels[index[j] - 1];

    }/*FOR*/
    /* compute the overall number of parameters as regressors plus intercept
     * and standard error (if continuous) or the number of levels - 1 (if
     * discrete) times the configurations of the discrete parents. */
    node_params = nconfig * (nlevels[i] == 0 ? ngp + 2 : nlevels[i] - 1);

    if (debugging)
      Rprintf("* node %s has %.0lf parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

    UNPROTECT(1);

  }/*FOR*/

  Free1D(nlevels);
  UNPROTECT(1);

  return ScalarReal(all_params);

}/*NPARAMS_CGNET*/


