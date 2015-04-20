#include "include/rcore.h"
#include "include/allocations.h"

/* get the number of parameters of the whole network (discrete case). */
SEXP nparams_dnet(SEXP graph, SEXP data, SEXP real, SEXP debug) {

int i = 0, j = 0, nnodes = 0;
int *index = NULL, r = isTRUE(real), debuglevel = isTRUE(debug);
double node_params = 0, all_params = 0;
double *nlevels = NULL;
SEXP nodes, node_data, parents, try;

  /* get nodes' number and data. */
  node_data = getListElement(graph, "nodes");
  nodes = getAttrib(node_data, R_NamesSymbol);
  nnodes = length(node_data);

  /* get the level count for each node. */
  nlevels = alloc1dreal(nnodes);
  for (i = 0; i < nnodes; i++)
    nlevels[i] = NLEVELS2(data, i);

  /* for each node... */
  for (i = 0; i < nnodes; i++) {

    /* reset the parameter counter. */
    node_params = 1;

    /* match the parents of the node. */
    parents = getListElement(VECTOR_ELT(node_data, i), "parents");
    PROTECT(try = match(nodes, parents, 0));
    index = INTEGER(try);

    /* compute the number of configurations. */
    for (j = 0; j < length(try); j++)
      node_params *= nlevels[index[j] - 1];

    UNPROTECT(1);

    /* multiply by the number of free parameters. */
    if (r > 0)
      node_params *= nlevels[i] - 1;
    else
      node_params *= nlevels[i];

    if (debuglevel > 0)
      Rprintf("* node %s has %.0lf parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

  }/*FOR*/

  return ScalarReal(all_params);

}/*NPARAMS_DNET*/

/* get the number of parameters of the whole network (continuous case). */
SEXP nparams_gnet(SEXP graph, SEXP debug) {

int i = 0, node_params = 0, debuglevel = isTRUE(debug);
double all_params = 0;
SEXP nodes = R_NilValue, temp = getListElement(graph, "nodes");

  if (debuglevel > 0)
    nodes = getAttrib(temp, R_NamesSymbol);

  /* add one parameter for each regressor, which means one for each
   * parent for each node plus the intercept. */
  for (i = 0; i < length(temp); i++) {

    node_params = length(getListElement(VECTOR_ELT(temp, i), "parents")) + 1;

    if (debuglevel > 0)
      Rprintf("* node %s has %d parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

  }/*FOR*/

  return ScalarReal(all_params);

}/*NPARAMS_GNET*/

/* get the number of parameters of the whole network (mixed case) */
SEXP nparams_cgnet(SEXP graph, SEXP data, SEXP debug) {

int i = 0, j = 0, nnodes = 0, debuglevel = isTRUE(debug);
int *nlevels = NULL, *index = NULL, ngp = 0;
double nconfig = 0, node_params = 0, all_params = 0;
SEXP nodes = R_NilValue, node_data, parents, temp;

  /* get nodes' number and data. */
  node_data = getListElement(graph, "nodes");
  nnodes = length(node_data);
  nodes = getAttrib(node_data, R_NamesSymbol);
  /* cache the number of levels of each variables (zero = continuous). */
  nlevels = Calloc(nnodes, int);
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
     * times configurations. */
    node_params = nconfig * (nlevels[i] == 0 ? ngp + 1 : nlevels[i]  - 1);

    if (debuglevel > 0)
      Rprintf("* node %s has %.0lf parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

    UNPROTECT(1);

  }/*FOR*/

  Free(nlevels);

  return ScalarReal(all_params);

}/*NPARAMS_CGNET*/

/* compute the number of parameters of a fitted model. */
SEXP fitted_nparams(SEXP bn, SEXP debug) {

int i = 0, j = 0, node_params = 0, nnodes = length(bn);
int res = 0, debuglevel = isTRUE(debug);
SEXP nodes = R_NilValue, node_data, temp;

  if (debuglevel > 0)
    nodes = getAttrib(bn, R_NamesSymbol);

  for (i = 0; i < nnodes; i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* get its probability distribution (if discrete). */
    temp = getListElement(node_data, "prob");

    if (!isNull(temp)) {

      /* reset the parameters' counter for this node. */
      node_params = 1;
      /* get the dimensions of the conditional probability table. */
      temp = getAttrib(temp, R_DimSymbol);
      /* compute the number of parameters. */
      for (j = 1; j < length(temp); j++)
        node_params *= INTEGER(temp)[j];

      node_params *= INTEGER(temp)[0] - 1;

    }/*THEN*/
    else {

      /* this is a continuous node, so it's a lot easier. */
      node_params = length(getListElement(node_data, "coefficients"));

    }/*ELSE*/

    if (debuglevel > 0)
      Rprintf("* node %s has %d parameter(s).\n", NODE(i), node_params);

    res += node_params;

  }/*FOR*/

  return ScalarInteger(res);

}/*FITTED_NPARAMS*/

