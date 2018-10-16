#include "include/rcore.h"
#include "include/sets.h"
#include "include/data.frame.h"
#include "include/scores.h"

/* posterior averagted Dirichlet probability (covers the BDLA score). */
double adpost(SEXP x, double l) {

int i = 0, a = 0, num = length(x);
int llx = NLEVELS(x), *xx = INTEGER(x), *n = NULL;
double imaginary = 0, alpha = 0, temp = 0, res = 0, iss_val = 0;

  /* initialize the contingency table. */
  n = Calloc1D(llx, sizeof(int));

  /* compute the frequency table of x, disregarding experimental data. */
  for (i = 0; i < num; i++)
    n[xx[i] - 1]++;

  for (a = 0; a < l; a++) {

    iss_val = R_pow(2, (1 - l) / 2 + a);

    imaginary = iss_val;
    alpha = imaginary / llx;

    /* compute the posterior probability. */
    for (i = 0, temp = 0; i < llx; i++)
      temp += lgammafn(n[i] + alpha) - lgammafn(alpha);
    temp += lgammafn(imaginary) - lgammafn(imaginary + num);

    res += temp / l;

  }/*FOR*/

  Free1D(n);

  return res;

}/*ADPOST*/

/* conditional averaged posterior Dirichlet probability (covers BDLA score). */
double acdpost(SEXP x, SEXP y, double l) {

int i = 0, j = 0, a = 0, num = length(x);
int llx = NLEVELS(x), lly = NLEVELS(y), p = llx * lly;
int *xx = INTEGER(x), *yy = INTEGER(y), **n = NULL, *nj = NULL;
double imaginary = 0, alpha = 0, temp = 0, temp2 = 0, res = 0, iss_val = 0;

  /* initialize the contingency table. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (i = 0; i < num; i++) {

    n[xx[i] - 1][yy[i] - 1]++;
    nj[yy[i] - 1]++;

  }/*FOR*/

  /* for every parents configuration ... */
  for (j = 0; j < lly; j++) {

    if (nj[j] == 0)
      continue;

    /* ... every imaginary sample size ... */
    for (a = 0, temp2 = 0; a < l; a++) {

      iss_val = R_pow(2, (1 - l) / 2 + a);

      imaginary = iss_val;
      alpha = imaginary / p;

      temp = lgammafn(imaginary / lly) - lgammafn(nj[j] + imaginary / lly);

      for (i = 0; i < llx; i++)
        temp += lgammafn(n[i][j] + alpha) - lgammafn(alpha);

      if (a == 0)
        temp2 = temp;
      else
        temp2 = logspace_add(temp2, temp);

    }/*FOR*/

    res += temp2 - log(l);

  }/*FOR*/

  Free1D(nj);
  Free2D(n, llx);

  return res;

}/*ACDPOST*/

/* Dirichlet averaged posterior probabilities (covers the BDLA score). */
double dirichlet_averaged_node(SEXP target, SEXP x, SEXP data, SEXP l,
    SEXP prior, SEXP beta, int sparse, int debuglevel) {

char *t = (char *)CHAR(STRING_ELT(target, 0));
double prob = 0, prior_prob = 0;
SEXP nodes, node_t, data_t, parents, parent_vars, config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
  /* compute the prior probability component for the node. */
  prior_prob = graph_prior_prob(prior, target, beta, nodes, debuglevel);

  if (length(parents) == 0) {

    prob = adpost(data_t, NUM(l));

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, !sparse));
    /* compute the marginal likelihood. */
    prob = acdpost(data_t, config, NUM(l));

    UNPROTECT(2);

  }/*ELSE*/

  if (debuglevel > 0) {

    Rprintf("  > (log)prior probability is %lf.\n", prior_prob);
    Rprintf("  > (log)posterior density is %lf.\n", prob);

  }/*THEN*/

  /* add the (log)prior to the marginal (log)likelihood to get the (log)posterior. */
  prob += prior_prob;

  UNPROTECT(1);

  return prob;

}/*DIRICHLET_AVERAGED_NODE*/

