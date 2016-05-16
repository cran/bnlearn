#include "include/rcore.h"
#include "include/sets.h"
#include "include/dataframe.h"
#include "include/scores.h"

/* posterior Dirichlet probability (covers BDe and K2 scores). */
double dpost(SEXP x, SEXP iss, SEXP exp) {

int i = 0, k = 0, num = length(x);
int llx = NLEVELS(x), *xx = INTEGER(x), *n = NULL, *imaginary = NULL;
double alpha = 0, res = 0;

  /* the correct vaules for the hyperparameters alpha are documented in
   * "Learning Bayesian Networks: The Combination of Knowledge and Statistical
   * Data" by Heckerman, Geiger & Chickering (1995), page 17. */

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distribution. */
    imaginary = &llx;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INTEGER(iss);
    alpha = (double) *imaginary / (double) llx;

  }/*ELSE*/

  /* initialize the contingency table. */
  n = Calloc1D(llx, sizeof(int));

  /* compute the frequency table of x, disregarding experimental data. */
  if (exp == R_NilValue) {

    for (i = 0; i < num; i++)
      n[xx[i] - 1]++;

  }/*THEN*/
  else {

    int *e = INTEGER(exp);

    for (i = 0, k = 0; i < num; i++) {
      if (i != e[k] - 1)
        n[xx[i] - 1]++;
      else
        k++;

    }/*FOR*/

    /* adjust the sample size to match the number of observational data. */
   num -= length(exp);

  }/*ELSE*/

  /* compute the posterior probability. */
  for (i = 0; i < llx; i++)
    res += lgammafn(n[i] + alpha) - lgammafn(alpha);
  res += lgammafn((double)(*imaginary)) -
            lgammafn((double)(*imaginary + num));

  Free1D(n);

  return res;

}/*DPOST*/

/* conditional posterior Dirichlet probability (covers BDe and K2 scores). */
double cdpost(SEXP x, SEXP y, SEXP iss, SEXP exp) {

int i = 0, j = 0, k = 0, imaginary = 0, num = length(x);
int llx = NLEVELS(x), lly = NLEVELS(y), p = llx * lly;
int *xx = INTEGER(x), *yy = INTEGER(y), **n = NULL, *nj = NULL;
double alpha = 0, res = 0;

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distribution. */
    imaginary = p;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INT(iss);
    alpha = ((double) imaginary) / ((double) p);

  }/*ELSE*/

  /* initialize the contingency table. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  if (exp == R_NilValue) {

    for (i = 0; i < num; i++) {

      n[xx[i] - 1][yy[i] - 1]++;
      nj[yy[i] - 1]++;

    }/*FOR*/

  }/*THEN*/
  else {

    int *e = INTEGER(exp);

    for (i = 0, k = 0; i < num; i++) {

      if (i != e[k] - 1) {

        n[xx[i] - 1][yy[i] - 1]++;
        nj[yy[i] - 1]++;

      }/*THEN*/
      else {

        k++;

      }/*ELSE*/

    }/*FOR*/

    /* adjust the sample size to match the number of observational data. */
   num -= length(exp);

  }/*ELSE*/

  /* compute the conditional posterior probability. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      res += lgammafn(n[i][j] + alpha) - lgammafn(alpha);
  for (j = 0; j < lly; j++)
    res += lgammafn((double)imaginary / lly) -
              lgammafn(nj[j] + (double)imaginary / lly);

  Free1D(nj);
  Free2D(n, llx);

  return res;

}/*CDPOST*/

/* Dirichlet posterior probabilities (covers BDe and K2 scores). */
double dirichlet_node(SEXP target, SEXP x, SEXP data, SEXP iss, SEXP prior,
    SEXP beta, SEXP experimental, int sparse, int debuglevel) {

char *t = (char *)CHAR(STRING_ELT(target, 0));
double prob = 0, prior_prob = 0;
SEXP nodes, node_t, data_t, exp_data, parents, parent_vars, config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  data_t = c_dataframe_column(data, target, TRUE, FALSE);
  /* extract the list of eperimental data. */
  exp_data = c_dataframe_column(experimental, target, TRUE, FALSE);
  /* compute the prior probability component for the node. */
  prior_prob = graph_prior_prob(prior, target, node_t, beta, debuglevel);

  if (length(parents) == 0) {

    prob = dpost(data_t, iss, exp_data);

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, !sparse));
    /* compute the marginal likelihood. */
    prob = cdpost(data_t, config, iss, exp_data);

    UNPROTECT(2);

  }/*ELSE*/

  if (debuglevel > 0) {

    Rprintf("  > (log)prior probability is %lf.\n", prior_prob);
    Rprintf("  > (log)posterior density is %lf.\n", prob);

  }/*THEN*/

  /* add the (log)prior to the marginal (log)likelihood to get the (log)posterior. */
  prob += prior_prob;

  return prob;

}/*DIRICHLET_NODE*/

