#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/dataframe.h"

double dlik(SEXP x, double *nparams) {

int i = 0;
int *n = NULL, *xx = INTEGER(x), llx = NLEVELS(x), num = length(x);
double res = 0;

  /* initialize the contingency table. */
  fill_1d_table(xx, &n, llx, num);

  /* compute the entropy from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    if (n[i] != 0)
      res += (double)n[i] * log((double)n[i] / num);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = llx - 1;

  Free1D(n);

  return res;

}/*DLIK*/

double cdlik(SEXP x, SEXP y, double *nparams) {

int i = 0, j = 0, k = 0;
int **n = NULL, *nj = NULL;
int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nj[j] += n[i][j];

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0)
        res += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (llx - 1) * lly;

  Free1D(nj);
  Free2D(n, llx);

  return res;

}/*CDLIK*/

double loglik_dnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel) {

double loglik = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, parent_vars, config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  data_t = c_dataframe_column(data, target, TRUE, FALSE);

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

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_DNODE*/
