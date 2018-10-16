#include "include/rcore.h"
#include "include/matrix.h"
#include "include/fitted.h"

/* get the number of parameters of the whole network (mixed case, also handles
 * discrete and Gaussian networks). */
SEXP nparams_cgnet(SEXP graph, SEXP data, SEXP debug) {

int i = 0, j = 0, nnodes = 0, debuglevel = isTRUE(debug);
int *nlevels = NULL, *index = NULL, ngp = 0;
double nconfig = 0, node_params = 0, all_params = 0;
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

    if (debuglevel > 0)
      Rprintf("* node %s has %.0lf parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

    UNPROTECT(1);

  }/*FOR*/

  Free1D(nlevels);
  UNPROTECT(1);

  return ScalarReal(all_params);

}/*NPARAMS_CGNET*/

/* compute the number of parameters of a fitted model. */
SEXP nparams_fitted(SEXP bn, SEXP effective, SEXP debug) {

int i = 0, j = 0, k = 0, node_params = 0, nnodes = length(bn), *pd = NULL;
int res = 0, debuglevel = isTRUE(debug), eff = isTRUE(effective);
double *ps = NULL, counter = 0;
SEXP nodes = R_NilValue, node_data, param_set, param_dims;
fitted_node_e node_type = ENOFIT;

  if (debuglevel > 0)
    PROTECT(nodes = getAttrib(bn, R_NamesSymbol));

  for (i = 0; i < nnodes; i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    node_type = r_fitted_node_label(node_data);

    switch (node_type) {

      case DNODE:
      case ONODE:
        /* get the conditional probability table and its dimension. */
        param_set = getListElement(node_data, "prob");
        param_dims = getAttrib(param_set, R_DimSymbol);
        pd = INTEGER(param_dims);
        ps = REAL(param_set);

        if (eff) {

          /* count the number of non-zero free parameters. */
          for (node_params = 0, k = 0; k < length(param_set) / pd[0]; k++) {

            /* check the elements of each conditional probability distribution. */
            for (counter = 0, j = 0; j < pd[0]; j++)
              counter += !ISNAN(ps[CMC(j, k, pd[0])]) && (ps[CMC(j, k, pd[0])] > 0);
            /* subtract the column total to get the free parameters. */
            if (counter > 0)
              counter--;

            node_params += counter;

          }/*FOR*/

        }/*THEN*/
        else {

          /* compute the number of parameters. */
          for (node_params = 1, j = 1; j < length(param_dims); j++)
            node_params *= pd[j];

          node_params *= pd[0] - 1;

        }/*ELSE*/

        break;

      case GNODE:
      case CGNODE:
        /* get the vector (or matrix) of regression coefficients. */
        param_set = getListElement(node_data, "coefficients");
        param_dims = getAttrib(param_set, R_DimSymbol);
        ps = REAL(param_set);

        if (eff) {

          /* count the number of non-zero regression coefficients. */
          for (node_params = 0, j = 0; j < length(param_set); j++)
            node_params += !ISNAN(ps[j]) && (ps[j] != 0);

        }/*THEN*/
        else {

          /* compute the number of parameters. */
          node_params = length(param_set);

        }/*ELSE*/

        /* residual standard errors are also parameters: Gaussian nodes have
         * just one, conditional Gaussian nodes have one for each configuration
         * of the discrete parents. */
        if (param_dims == R_NilValue)
          node_params += 1;
        else
          node_params += INTEGER(param_dims)[1];

        break;

      default:
        error("unknown node type (class: %s).",
           CHAR(STRING_ELT(getAttrib(node_data, R_ClassSymbol), 0)));

    }/*SWITCH*/

    if (debuglevel > 0)
      Rprintf("* node %s has %d parameter(s).\n", NODE(i), node_params);

    res += node_params;

  }/*FOR*/

  if (debuglevel > 0)
    UNPROTECT(1);

  return ScalarInteger(res);

}/*NPARAMS_FITTED*/

