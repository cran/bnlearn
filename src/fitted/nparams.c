#include "../include/rcore.h"
#include "../fitted/fitted.h"
#include "../math/linear.algebra.h"

/* compute the number of parameters of a fitted model. */
double nparams_fitted_node(ldist ld, fitted_node_e type, bool effective) {

int j = 0, k = 0, nrows = 0, nconfigs = 0, counter = 0;
double node_params = 0;

  switch(type) {

    /* ... if it's a discrete node: */
    case DNODE:
    case ONODE:

      nrows = ld.d.dims[0];
      nconfigs = ld.d.nconfigs;

      if (effective) {

        for (j = 0; j < nconfigs; j++) {

          /* ... and count how many conditional probabilities:
           * 1) are not NA;
           * 2) are not zeroes, as suggested by Fienberg (controversial!); */
          for (k = 0, counter = 0; k < nrows; k++)
            counter += !ISNAN(ld.d.cpt[CMC(k, j, nrows)]) &&
                       (ld.d.cpt[CMC(k, j, nrows)] > 0);
          /* 3) are not bound by the column totals. */
          if (counter > 0)
            counter--;

          /* singular conditional distributions have zero parameters. */

          node_params += counter;

        }/*FOR*/

      }/*THEN*/
      else {

        /* the number of parameters for each node is (number of rows - 1),
         * either in total or for each parents configuration. */
        node_params = nrows * nconfigs - nconfigs;

      }/*ELSE*/

      break;

    /* ... if it's a Gaussian node: */
    case GNODE:

      if (effective) {

        /* count the number of non-zero regression coefficients... */
        for (j = 0; j < ld.g.ncoefs; j++)
          node_params += !ISNAN(ld.g.coefs[j]) && (ld.g.coefs[j] != 0);
        /* ... and the standard error of the residuals. */
        node_params++;

      }/*THEN*/
      else {

        /*... the number of parameters is (number of regression coefficients,
         * one per parents) plus two (intercept and standard error). */
        node_params = ld.nparents + 2;

      }/*ELSE*/

      break;

    /* ... if it's a conditional Gaussian node: */
    case CGNODE:

      nconfigs = ld.cg.nconfigs;
      nrows = ld.cg.ncoefs;

      if (effective) {

        /* count the number of non-zero regression coefficients... */
        for (j = 0; j < nconfigs * nrows; j++)
          node_params += !ISNAN(ld.cg.coefs[j]) && (ld.cg.coefs[j] != 0);
        /* ... and all the standard errors of the residuals (one for each
         * configuration of the discrete parents). */
        node_params += nconfigs;

      }/*THEN*/
      else {

        /* the number of parameters is (number of continuous parents + 2) for
         * each configuration of the discrete parents (of which there should
         * be at least one). */
        node_params += nconfigs * (nrows + 1);

      }/*ELSE*/

      break;

    case ENOFIT:
    default:
      break;

  }/*SWITCH*/

  return node_params;

}/*NPARAMS_FITTED_NODE*/

