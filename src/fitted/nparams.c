#include "../include/rcore.h"
#include "../fitted/fitted.h"
#include "../math/linear.algebra.h"

/* compute the number of parameters of a fitted model. */
double nparams_fitted_node(ldist ld, fitted_node_e type) {

int nrows = 0, nconfigs = 0;
double node_params = 0;

  switch(type) {

    /* ... if it's a discrete node: */
    case DNODE:
    case ONODE:

      nrows = ld.d.dims[0];
      nconfigs = ld.d.nconfigs;

      /* the number of parameters for each node is (number of rows - 1),
       * either in total or for each parents configuration. */
      node_params = nrows * nconfigs - nconfigs;

      break;

    /* ... if it's a Gaussian node: */
    case GNODE:

      /* the number of parameters is (number of regression coefficients, one
       * per parent) plus two (intercept and standard error). */
      node_params = ld.nparents + 2;

      break;

    /* ... if it's a conditional Gaussian node: */
    case CGNODE:

      nconfigs = ld.cg.nconfigs;
      nrows = ld.cg.ncoefs;

      /* the number of parameters is (number of continuous parents + 2) for each
       * configuration of the discrete parents (of which there should be at
       * least one). */
      node_params += nconfigs * (nrows + 1);

      break;

    case ENOFIT:
    default:
      break;

  }/*SWITCH*/

  return node_params;

}/*NPARAMS_FITTED_NODE*/

