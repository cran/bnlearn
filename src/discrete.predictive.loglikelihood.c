#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/data.frame.h"
#include "include/data.table.h"

double pdnode(SEXP x, SEXP new_x, double *nparams) {

int i = 0, num = length(x), num2 = length(new_x);
int *n = NULL, *n2 = NULL, *xx = INTEGER(x), *xx2 = INTEGER(new_x), llx = NLEVELS(x);
double res = 0;

  /* initialize the contingency table. */
  fill_1d_table(xx, &n, llx, num);
  fill_1d_table(xx2, &n2, llx, num2);

  /* compute the entropy from the joint and marginal frequencies. */
  for (i = 0; i < llx; i++)
    if (n[i] != 0)
      res += (double)n2[i] * log((double)n[i] / num2);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = llx - 1;

  Free1D(n);
  Free1D(n2);

  return res;

}/*PDNODE*/

double cpdnode(SEXP x, SEXP y, SEXP x2, SEXP y2, double *nparams) {

int i = 0, j = 0, k = 0;
int **n = NULL, **n2 = NULL, *nj = NULL;
int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x), num2 = length(x2);
int *xx = INTEGER(x), *yy = INTEGER(y), *xx2 = INTEGER(x2), *yy2 = INTEGER(y2);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  n2 = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;
  for (k = 0; k < num2; k++)
    n2[xx2[k] - 1][yy2[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nj[j] += n[i][j];

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0)
        res += (double)n2[i][j] * log((double)n[i][j] / (double)nj[j]);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (llx - 1) * lly;

  Free1D(nj);
  Free2D(n, llx);
  Free2D(n2, llx);

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
  data_t = c_dataframe_column(data, target, TRUE, FALSE);
  new_t = c_dataframe_column(newdata, target, TRUE, FALSE);

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

  return loglik;

}/*PREDICTIVE_LOGLIK_DNODE*/

