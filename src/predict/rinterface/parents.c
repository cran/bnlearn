#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/data.table.h"
#include "../../core/sets.h"
#include "../../fitted/fitted.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/data.frame.h"
#include "../predict.h"

/* predict the values of a node from its parents, by node type. */
SEXP predict_from_parents(SEXP fitted, SEXP node, SEXP data, SEXP prob,
    SEXP debug) {

bool debugging = isTRUE(debug), include_prob = isTRUE(prob);
fitted_bn bn = fitted_network_from_SEXP(fitted);
const char *node_name = CHAR(STRING_ELT(node, 0));
int node_id = 0, nobs = 0;
fitted_node_e type;
ldist ld;
tabular par_tab = { 0 };
SEXP temp, par_data, result, probtab = R_NilValue;

  /* locate the target node in the fitted network. */
  for (node_id = 0; node_id < bn.nnodes; node_id++)
    if (strcmp(bn.labels[node_id], node_name) == 0)
      break;
  type = bn.node_types[node_id];
  ld = bn.ldists[node_id];

  /* extract the parents of the node from the fitted network... */
  temp = getListElement(fitted, (char *)node_name);
  temp = getListElement(temp, "parents");
  /* ... subset them out of the data (this aligns the columns with the
   * parameters) and transform them into a tabular. */
  PROTECT(par_data = c_dataframe_column(data, temp, FALSE, FALSE));
  par_tab = tabular_from_SEXP(par_data, 0, 0);
  nobs = par_tab.m.nobs;

  /* allocate the return value. */
  PROTECT(result = fitnode2df(fitted, STRING_ELT(node, 0), nobs));

  /* dispatch by node type. */
  switch(type) {

    case GNODE:

      c_gaussian_predict(ld, par_tab.ccol, par_tab.nccols, nobs, REAL(result),
        debugging);

      break;

    case DNODE:
    case ONODE: {

      int *configs = NULL, ncfg = 0;
      double *pt = NULL;

      /* allocate the table of the prediction probabilities. */
      if (include_prob) {

        PROTECT(probtab = allocMatrix(REALSXP, ld.d.dims[0], nobs));
        pt = REAL(probtab);
        memset(pt, '\0', nobs * ld.d.dims[0] * sizeof(double));

      }/*THEN*/

      /* collapse the discrete parents into a configuration index. */
      if (par_tab.ndcols > 0) {

        configs = Calloc1D(nobs, sizeof(int));
        c_fast_config(par_tab.dcol, nobs, par_tab.ndcols, par_tab.nlvl, configs,
          &ncfg, 1);

      }/*THEN*/

      c_discrete_predict(ld, configs, ld.d.nconfigs, nobs, INTEGER(result), pt,
        include_prob, debugging);

      Free1D(configs);

      break;

    }/*DNODE*/

    case CGNODE: {

      int *configs = Calloc1D(nobs, sizeof(int)), ncfg = 0;

      /* collapse the discrete parents into a configuration index (a single
       * configuration shared by all the observations if there are none). */
      if (par_tab.ndcols > 0)
        c_fast_config(par_tab.dcol, nobs, par_tab.ndcols, par_tab.nlvl, configs,
          &ncfg, 1);
      else
        for (int i = 0; i < nobs; i++)
          configs[i] = 1;

      c_cgaussian_predict(ld, par_tab.ccol, par_tab.nccols, configs, nobs,
        REAL(result), debugging);

      Free1D(configs);

      break;

    }/*CGNODE*/

    case ZIHPNODE:
    case ZINBNODE:

      c_zero_inflated_predict(ld, type, par_tab, nobs, REAL(result), debugging);

      break;

    default:
      FreeTAB(par_tab);
      FreeFittedBN(bn);
      error("unknown node type.");

  }/*SWITCH*/

  /* attach the prediction probabilities to discrete predictions. */
  if (probtab != R_NilValue) {

    setDimNames(probtab, getAttrib(result, R_LevelsSymbol), R_NilValue);
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(3);

  }/*THEN*/
  else {

    UNPROTECT(2);

  }/*ELSE*/

  FreeTAB(par_tab);
  FreeFittedBN(bn);

  return result;

}/*PREDICT_FROM_PARENTS*/
