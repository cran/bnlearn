#include "common.h"
#include <Rmath.h>

double c_gloss(int *cur, SEXP cur_parents, double *coefs, double *sd,
    double **columns, SEXP nodes, int ndata, double *per_sample);
double c_dloss(int *cur, SEXP cur_parents, int *configs, double *prob,
    SEXP data, SEXP nodes, int ndata, int nlevels, double *per_sample);

#define DISCRETE 0
#define GAUSSIAN 1

SEXP entropy_loss(SEXP fitted, SEXP orig_data, SEXP by_sample, SEXP keep,
    SEXP debug) {

int i = 0, k = 0, ndata = 0, nnodes = length(fitted), nlevels = 0, type = 0;
int *configs = NULL, *debuglevel = LOGICAL(debug), *by = LOGICAL(by_sample);
int *to_keep = NULL;
double *res = 0, *res_sample = NULL, **columns = 0, cur_loss = 0;
const char *class = NULL;
SEXP data, cur_node, nodes, result, result_sample, coefs, sd, parents, try;

  /* get the node labels. */
  nodes = getAttrib(fitted, R_NamesSymbol);
  /* rearrange the columns of the data to match the network. */
  PROTECT(data = c_dataframe_column(orig_data, nodes, FALSE, TRUE));
  /* get the sample size. */
  ndata = length(VECTOR_ELT(data, 0));
  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;
  /* allocate the sample's contributions if needed. */
  if (*by > 0) {

    PROTECT(result_sample = allocVector(REALSXP, ndata));
    res_sample = REAL(result_sample);
    memset(res_sample, '\0', ndata * sizeof(double));

  }/*THEN*/

  /* find out which nodes to use in computing the entropy loss. */
  PROTECT(try = match(nodes, keep, 0));
  to_keep = INTEGER(try);
  R_isort(to_keep, length(try));

  /* determine the class of the fitted network. */
  class = CHAR(STRING_ELT(getAttrib(VECTOR_ELT(fitted, 0), R_ClassSymbol), 0));

  if (strcmp(class, "bn.fit.gnode") == 0) {

    /* dereference the data set's columns. */
    columns = (double **) alloc1dpointer(nnodes);
    for (i = 0; i < nnodes; i++)
      columns[i] = REAL(VECTOR_ELT(data, i));

    type = GAUSSIAN;

  }/*THEN*/
  else if ((strcmp(class, "bn.fit.dnode") == 0) || (strcmp(class, "bn.fit.onode") == 0)) {

    /* allocate an array for parents' configurations. */
    configs = alloc1dcont(ndata);

    type = DISCRETE;

  }/*THEN*/

  /* iterate over the nodes. */
  for (i = 0; i < nnodes; i++) {

    if (i == to_keep[k] - 1) {

      k++;

    }/*THEN*/
    else {

      if (*debuglevel > 0)
        Rprintf("  > skipping node %s.\n", NODE(i));

      continue;

    }/*ELSE*/

    /* get the current node. */
    cur_node = VECTOR_ELT(fitted, i);
    /* get the parents of the node. */
    parents = getListElement(cur_node, "parents");
    /* get the parameters (regression coefficients and residuals' standard
     * deviation for Gaussian nodes, conditional probabilities for discrete
     * nodes), and compute the loss. */
    switch(type)  {

      case GAUSSIAN:

        coefs = getListElement(cur_node, "coefficients");
        sd = getListElement(cur_node, "sd");

        cur_loss = c_gloss(&i, parents, REAL(coefs), REAL(sd), columns, nodes,
                     ndata, res_sample);
        break;

      case DISCRETE:

        coefs = getListElement(cur_node, "prob");
        nlevels = INT(getAttrib(coefs, R_DimSymbol));

        cur_loss = c_dloss(&i, parents, configs, REAL(coefs), data, nodes,
                     ndata, nlevels, res_sample);
        break;

    }/*SWITCH*/

    if (*debuglevel > 0)
      Rprintf("  > log-likelihood loss for node %s is %lf.\n", NODE(i), cur_loss);

    /* add the node contribution to the return value. */
    *res += cur_loss;

  }/*FOR*/

  if (*by > 0) {

    UNPROTECT(4);
    return result_sample;

  }/*THEN*/
  else {

    UNPROTECT(3);
    return result;

  }/*ELSE*/

}/*ENTROPY_LOSS*/

/* Gaussian loss for a single node. */
double c_gloss(int *cur, SEXP cur_parents, double *coefs, double *sd,
    double **columns, SEXP nodes, int ndata, double *per_sample) {

int i = 0, j = 0, *p = NULL, nparents = length(cur_parents);
double mean = 0, logprob = 0, result = 0;
SEXP try;

  if (nparents > 0) {

    PROTECT(try = match(nodes, cur_parents, 0));
    p = INTEGER(try);

  }/*THEN*/

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    mean = coefs[0];

    for (j = 0; j < nparents; j++)
      mean += columns[p[j] - 1][i] * coefs[j + 1];

    /* compute the log-likelihood of this observation. */
    logprob = dnorm(columns[*cur][i], mean, *sd, TRUE);

    result += logprob;

    if (per_sample)
      per_sample[i] += logprob;

  }/*FOR*/

  if (nparents > 0)
    UNPROTECT(1);

  /* switch to the negentropy. */
  result /= -ndata;

  return result;

}/*C_GLOSS*/

/* multinomial loss for a single node. */
double c_dloss(int *cur, SEXP cur_parents, int *configs, double *prob,
    SEXP data, SEXP nodes, int ndata, int nlevels, double *per_sample) {

int i = 0, dropped = 0, *obs = NULL;
double logprob = 0, result = 0;
SEXP temp_df;

  /* get the target variable. */
  obs = INTEGER(VECTOR_ELT(data, *cur));
  /* get the parents' configurations. */
  if (length(cur_parents) > 0) {

    PROTECT(temp_df = c_dataframe_column(data, cur_parents, FALSE, FALSE));
    cfg(temp_df, configs, NULL);

    for (i = 0; i < ndata; i++) {

      logprob = log(prob[CMC(obs[i] - 1, configs[i], nlevels)]);

      if (!R_FINITE(logprob) || ISNAN(logprob))
        dropped++;
      else
        result += logprob;

      if (per_sample)
        per_sample[i] += logprob;

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else {

    for (i = 0; i < ndata; i++) {

      logprob = log(prob[obs[i] - 1]);

      if (!R_FINITE(logprob) || ISNAN(logprob))
        dropped++;
      else
        result += logprob;

      if (per_sample)
        per_sample[i] += logprob;

    }/*FOR*/

  }/*ELSE*/

  /* switch to the negentropy. */
  result /= -(ndata - dropped);

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because the corresponding probabilities for node %s were 0 or NaN.", dropped, NODE(*cur));

  return result;

}/*C_DLOSS*/

/* classification error of a single node as a loss function. */
SEXP class_err(SEXP reference, SEXP predicted) {

int i = 0, dropped = 0, ndata = length(reference);
int *r = INTEGER(reference), *p = INTEGER(predicted);
double *res = NULL;
SEXP result;

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* count how many elements differ (this assumes the levels of the two factors
   * are the same and in the same order); NAs are dropped. */
  for (i = 0; i < ndata; i++) {

    if ((r[i] == NA_INTEGER) || (p[i] == NA_INTEGER))
      dropped++;
    else if ((r[i] != p[i]))
      (*res)++;

  }/*FOR*/

  /* rescale into a probability. */
  *res /= (ndata - dropped);

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because of missing values.", dropped);

  UNPROTECT(1);

  return result;

}/*CLASS_ERR*/
