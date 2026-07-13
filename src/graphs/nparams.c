#include "graphs.h"

/* number of parameters of a single node. */
double nparams_node_count(sparse_amat parents, tabular dt, int i,
    estimator_e estimator) {

int ncp = 0;
double nconfig = 1, node_params = 0;

  /* count the continuous parents (regressors) and the configurations of the
   * discrete parents, reading the parent indexes straight from the CSR row. */
  for (int k = parents.rowptr[i]; k < parents.rowptr[i + 1]; k++) {

    int p = parents.colidx[k];
    if (dt.m.flag[p].continuous)
      ncp++;
    else
      nconfig *= dt.nlvl[dt.map[p]];

  }/*FOR*/

  if (dt.m.flag[i].discrete) {

    /* discrete node: (levels - 1) free parameters. */
    node_params = dt.nlvl[dt.map[i]] - 1;

  }/*THEN*/
  else {

    /* continuous node: the count depends on the estimator. */
    switch (estimator) {

      /* zero-inflated count nodes: (2 * regressors) + 2 intercepts + dispersion. */
      case MLE_ZIHP:
      case MLE_ZINB:
        node_params = 2 * ncp + 3;
        break;

      /* Gaussian nodes: regressors + intercept + standard error. */
      case MLE_G:
      case HARD_EM_G:
      case MLE_CG:
      case HARD_EM_CG:
      default:
        node_params = ncp + 2;
        break;

    }/*SWITCH*/

  }/*ELSE*/

  /* the count is replicated once per configuration of the discrete parents. */
  return nconfig * node_params;

}/*NPARAMS_NODE_COUNT*/
