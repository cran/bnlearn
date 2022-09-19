#include "../include/rcore.h"
#include "../core/sets.h"
#include "../minimal/data.frame.h"
#include "../minimal/common.h"
#include "../core/contingency.tables.h"
#include "scores.h"

static double cdlik_diff(counts2d joint) {

double res = 0;

  for (int j = 0; j < joint.lly; j++) {

    for (int i = 0; i < joint.llx; i++)
      if (joint.n[i][j] != 0)
        res += (double)joint.n[i][j] * log((double)joint.n[i][j]);

    if (joint.nj[j] != 0)
      res -= (double)joint.nj[j] * log((double)joint.nj[j]);

  }/*FOR*/

  return res;

}/*CDLIK_DIFF*/

/* factorized normalized maximum likelihood, root nodes. */
double fnml(SEXP x) {

  double res = 0;
  counts1d marginal = { 0 };

  /* initialize the contingency table. */
  marginal = new_1d_table(NLEVELS(x));
  fill_1d_table(INTEGER(x), &marginal, length(x));
  res = dlik(marginal);
  res -= nml_regret(marginal.nobs, NLEVELS(x));
  Free1DTAB(marginal);

  return res;

}/*FNML*/

/* factorized normalized maximum likelihood, nodes with parents. */
double cfnml(SEXP x, SEXP y) {

  double res = 0;
  counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(NLEVELS(x), NLEVELS(y), TRUE);
  fill_2d_table(INTEGER(x), INTEGER(y), &joint, length(x));

  res = cdlik(joint);
  for (int j = 0; j < joint.lly; j++) {
    if (joint.nj[j] != 0)
      res -= nml_regret(joint.nj[j], NLEVELS(x));
  }
  Free2DTAB(joint);

  return res;

}/*CFNML*/

double fnml_node(SEXP target, SEXP x, SEXP data, bool debugging) {

double nml = 0;
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

    nml = fnml(data_t);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    nml = cfnml(data_t, config);

    UNPROTECT(2);

  }/*ELSE*/

  if (debugging)
    Rprintf("  > normalized maximum likelihood is %lf.\n", nml);

  UNPROTECT(1);

  return nml;

}/*FNML_NODE*/


/* quotient normalized maximum likelihood, root nodes. */
double qnml(SEXP x) {

  return fnml(x);

}/*QNML*/

/* quotient normalized maximum likelihood, nodes with parents. */
double cqnml(SEXP x, SEXP y) {

  double res = 0;
  counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(NLEVELS(x), NLEVELS(y), TRUE);
  fill_2d_table(INTEGER(x), INTEGER(y), &joint, length(x));
  res += cdlik_diff(joint);
  res -= nml_regret(joint.nobs, joint.llx * joint.lly);
  res += nml_regret(joint.nobs, joint.lly);
  Free2DTAB(joint);

  return res;

}/*CQNML*/

double qnml_node(SEXP target, SEXP x, SEXP data, bool debugging) {

double nml = 0;
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

    nml = qnml(data_t);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    nml = cqnml(data_t, config);

    UNPROTECT(2);

  }/*ELSE*/

  if (debugging)
    Rprintf("  > normalized maximum likelihood is %lf.\n", nml);

  UNPROTECT(1);

  return nml;

}/*QNML_NODE*/
