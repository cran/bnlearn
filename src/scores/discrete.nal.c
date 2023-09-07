#include "../include/rcore.h"
#include "../core/sets.h"
#include "../minimal/data.frame.h"
#include "../minimal/common.h"
#include "../core/contingency.tables.h"
#include "scores.h"

double nal_dnode_root(SEXP x, double k) {

double res = 0;
counts1d marginal = { 0 };

  /* initialize the contingency table. */
  marginal = new_1d_table(NLEVELS(x));
  fill_1d_table(INTEGER(x), &marginal, length(x));

  if (marginal.nobs == 0) {

    /* no locally-complete observations: prevent a network with this local
     * distribution from being selected in structure learning. */
    res = R_NegInf;

  }/*THEN*/
  else {

    /* compute the node-average log-likelihood from the marginal frequencies. */
    res = dlik(marginal) / marginal.nobs;
    /* add the penalty term, scaled by the original sample size. */
    res -= k / length(x) * (marginal.llx - 1);

  }/*ELSE*/

  Free1DTAB(marginal);

  return res;

}/*NAL_DNODE_ROOT*/

double nal_dnode_parents(SEXP x, SEXP y, double k) {

double res = 0;
counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(NLEVELS(x), NLEVELS(y), TRUE);
  fill_2d_table(INTEGER(x), INTEGER(y), &joint, length(x));

  if (joint.nobs == 0) {

    /* no locally-complete observations: prevent a network with this local
     * distribution from being selected in structure learning. */
    res = R_NegInf;

  }/*THEN*/
  else {

    /* compute the node-average log-likelihood from the marginal frequencies. */
    res = cdlik(joint) / joint.nobs;
    /* add the penalty term, scaled by the original sample size. */
    res -= k / length(x) * ((joint.llx - 1) * joint.lly);

  }/*ELSE*/

  Free2DTAB(joint);

  return res;

}/*NAL_DNODE_PARENTS*/

double nal_dnode(SEXP target, SEXP x, SEXP data, double k, bool debugging) {

double loglik = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, parent_vars, config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));

  if (length(parents) == 0) {

    loglik = nal_dnode_root(data_t, k);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    loglik = nal_dnode_parents(data_t, config, k);

    UNPROTECT(2);

  }/*ELSE*/

  if (debugging)
    Rprintf("  > log-likelihood is %lf.\n", loglik);

  UNPROTECT(1);

  return loglik;

}/*NAL_DNODE*/
