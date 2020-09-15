#include "include/rcore.h"
#include "include/sets.h"
#include "include/contingency.tables.h"
#include "include/data.frame.h"
#include "include/data.table.h"

double pdnode(SEXP x, SEXP new_x, double *nparams) {

double res = 0;
counts1d train = { 0 }, test = { 0 };

  /* initialize the contingency tables for train and test data. */
  train = new_1d_table(NLEVELS(x));
  test = new_1d_table(train.llx);
  fill_1d_table(INTEGER(x), &train, length(x));
  fill_1d_table(INTEGER(new_x), &test, length(new_x));

  /* compute the entropy from the joint and marginal frequencies. */
  for (int i = 0; i < train.llx; i++)
    if (train.n[i] != 0)
      res += (double)test.n[i] * log((double)train.n[i] / test.nobs);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = train.llx - 1;

  Free1DTAB(train);
  Free1DTAB(test);

  return res;

}/*PDNODE*/

double cpdnode(SEXP x, SEXP y, SEXP x2, SEXP y2, double *nparams) {

double res = 0;
counts2d train = { 0 }, test = { 0 };

  /* initialize the contingency tables for train and test data. */
  train = new_2d_table(NLEVELS(x), NLEVELS(y), TRUE);
  test = new_2d_table(train.llx, train.lly, FALSE);
  fill_2d_table(INTEGER(x), INTEGER(y), &train, length(x));
  fill_2d_table(INTEGER(x2), INTEGER(y2), &test, length(x2));

  /* compute the conditional entropy from the joint and marginal frequencies. */
  for (int i = 0; i < train.llx; i++)
    for (int j = 0; j < train.lly; j++)
      if (train.n[i][j] != 0)
        res += (double)test.n[i][j] * log((double)train.n[i][j] / train.nj[j]);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (train.llx - 1) * train.lly;

  Free2DTAB(train);
  Free2DTAB(test);

  return res;

}/*CPDNODE*/

double predictive_loglik_dnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, int debuglevel) {

double loglik = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, new_t;
SEXP parent_vars, new_parents, config, new_config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
  PROTECT(new_t = c_dataframe_column(newdata, target, TRUE, FALSE));

  if (length(parents) == 0) {

    loglik = pdnode(data_t, new_t, nparams);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    PROTECT(new_parents = c_dataframe_column(newdata, parents, FALSE, FALSE));
    PROTECT(new_config = c_configurations(new_parents, TRUE, TRUE));
    /* compute the log-likelihood. */
    loglik = cpdnode(data_t, config, new_t, new_config, nparams);

    UNPROTECT(4);

  }/*ELSE*/

  //loglik -= /* log(test_size)/2 */ (*nparams);

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  UNPROTECT(2);

  return loglik;

}/*PREDICTIVE_LOGLIK_DNODE*/

