#include "include/rcore.h"
#include "include/sets.h"
#include "include/data.frame.h"
#include "include/contingency.tables.h"

double dlik(SEXP x, double *nparams) {

double res = 0;
counts1d marginal = { 0 };

  /* initialize the contingency table. */
  marginal = new_1d_table(NLEVELS(x));
  fill_1d_table(INTEGER(x), &marginal, length(x));

  /* compute the entropy from the marginal frequencies. */
  for (int i = 0; i < marginal.llx; i++)
    if (marginal.n[i] != 0)
      res += (double)marginal.n[i] * log((double)marginal.n[i] / marginal.nobs);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = marginal.llx - 1;

  Free1DTAB(marginal);

  return res;

}/*DLIK*/

double cdlik(SEXP x, SEXP y, double *nparams) {

double res = 0;
counts2d joint = { 0 };

  /* initialize the contingency table and the marginal frequencies. */
  joint = new_2d_table(NLEVELS(x), NLEVELS(y), TRUE);
  fill_2d_table(INTEGER(x), INTEGER(y), &joint, length(x));

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (int i = 0; i < joint.llx; i++)
    for (int j = 0; j < joint.lly; j++)
      if (joint.n[i][j] != 0)
        res += (double)joint.n[i][j] * log((double)joint.n[i][j] / joint.nj[j]);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (joint.llx - 1) * joint.lly;

  Free2DTAB(joint);

  return res;

}/*CDLIK*/

double loglik_dnode(SEXP target, SEXP x, SEXP data, double *nparams,
    bool debugging) {

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

    loglik = dlik(data_t, nparams);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    loglik = cdlik(data_t, config, nparams);

    UNPROTECT(2);

  }/*ELSE*/

  if (debugging)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  UNPROTECT(1);

  return loglik;

}/*LOGLIK_DNODE*/
