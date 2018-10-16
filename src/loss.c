#include "include/rcore.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/sets.h"
#include "include/fitted.h"

double c_gloss(int *cur, SEXP cur_parents, double *coefs, double *sd,
    void **columns, SEXP nodes, int ndata, double *per_sample,
    int allow_singular);
double c_dloss(int *cur, SEXP cur_parents, int *configs, double *prob,
    SEXP data, SEXP nodes, int ndata, int nlevels, double *per_sample,
    int *dropped);
double c_cgloss(int *cur, SEXP cur_parents, SEXP dparents, SEXP gparents,
    SEXP dlevels, double *coefs, double *sd, void **columns, SEXP nodes,
    int ndata, double *per_sample, int allow_singular, int *dropped);
double c_entropy_loss(SEXP fitted, SEXP orig_data, int ndata, int by,
    double *res_sample, SEXP keep, int allow_singular, int warnlevel,
    int debuglevel);

#define UPDATE_LOSS(cond) \
      if (cond) \
        (*dropped)++; \
      else \
        result += logprob; \
      if (per_sample) \
        per_sample[i] += logprob;

SEXP entropy_loss(SEXP fitted, SEXP data, SEXP by_sample, SEXP keep,
    SEXP debug) {

int *by = LOGICAL(by_sample), ndata = length(VECTOR_ELT(data, 0));
double *res_sample = NULL, loss = 0;
SEXP result_sample = R_NilValue;

  /* allocate the sample's contributions if needed. */
  if (*by) {

    PROTECT(result_sample = allocVector(REALSXP, ndata));
    res_sample = REAL(result_sample);
    memset(res_sample, '\0', ndata * sizeof(double));

  }/*THEN*/

  loss = c_entropy_loss(fitted, data, ndata, *by, res_sample, keep,
                  TRUE, TRUE, isTRUE(debug));

  if (*by)
    UNPROTECT(1);

  return (*by) ? result_sample : ScalarReal(loss);

}/*ENTROPY_LOSS*/

double c_entropy_loss(SEXP fitted, SEXP orig_data, int ndata, int by,
    double *res_sample, SEXP keep, int allow_singular, int warnlevel,
    int debuglevel) {

int i = 0, k = 0, nnodes = length(fitted), nlevels = 0, dropped = 0;
int *configs = NULL, *to_keep = NULL;
double result = 0, cur_loss = 0;
fitted_node_e node_type = ENOFIT;
void **columns = NULL;
SEXP data, cur_node, nodes, coefs, sd, parents, try;
SEXP dparents, gparents, dlevels;

  /* get the node labels. */
  PROTECT(nodes = getAttrib(fitted, R_NamesSymbol));
  /* rearrange the columns of the data to match the network. */
  PROTECT(data = c_dataframe_column(orig_data, nodes, FALSE, TRUE));
  /* find out which nodes to use in computing the entropy loss. */
  PROTECT(try = match(nodes, keep, 0));
  to_keep = INTEGER(try);
  R_isort(to_keep, length(try));

  /* dereference the data set's columns. */
  columns = Calloc1D(nnodes, sizeof(void *));
  for (i = 0; i < nnodes; i++)
    columns[i] = (void *) DATAPTR(VECTOR_ELT(data, i));

  /* allocate an array for parents' configurations. */
  if (!c_is(fitted, "bn.fit.gnet"))
    configs = Calloc1D(ndata, sizeof(int));

  /* iterate over the nodes. */
  for (i = 0, k = 0; i < nnodes; i++) {

    if (i == to_keep[k] - 1) {

      /* prevent k from overflowing but do not break out of the loop, to allow
       * the debuggin output to cover all the nodes. */
      if (k < length(try) - 1)
        k++;

    }/*THEN*/
    else {

      if (debuglevel > 0)
        Rprintf("  > skipping node %s.\n", NODE(i));

      continue;

    }/*ELSE*/

    /* reset the counter for the dropped observations. */
    dropped = 0;
    /* get the current node. */
    cur_node = VECTOR_ELT(fitted, i);
    /* get the parents of the node. */
    parents = getListElement(cur_node, "parents");
    /* get the parameters (regression coefficients and residuals' standard
     * deviation for Gaussian nodes, conditional probabilities for discrete
     * nodes), and compute the loss. */
    node_type = r_fitted_node_label(cur_node);

    switch(node_type) {

      case DNODE:
      case ONODE:
        coefs = getListElement(cur_node, "prob");
        nlevels = INT(getAttrib(coefs, R_DimSymbol));

        cur_loss = c_dloss(&i, parents, configs, REAL(coefs), data, nodes,
                     ndata, nlevels, res_sample, &dropped);

        break;


      case GNODE:
        coefs = getListElement(cur_node, "coefficients");
        sd = getListElement(cur_node, "sd");

        cur_loss = c_gloss(&i, parents, REAL(coefs), REAL(sd), columns, nodes,
                     ndata, res_sample, allow_singular);

        break;

      case CGNODE:

        coefs = getListElement(cur_node, "coefficients");
        sd = getListElement(cur_node, "sd");
        dparents = getListElement(cur_node, "dparents");
        gparents = getListElement(cur_node, "gparents");
        dlevels = getListElement(cur_node, "dlevels");

        cur_loss = c_cgloss(&i, parents, dparents, gparents, dlevels,
                     REAL(coefs), REAL(sd), columns, nodes, ndata, res_sample,
                     allow_singular, &dropped);

        break;

      default:
        error("unknown node type (class: %s).",
           CHAR(STRING_ELT(getAttrib(cur_node, R_ClassSymbol), 0)));

    }/*SWITCH*/

    /* print a warning if data were dropped. */
    if ((warnlevel > 0) && (dropped > 0))
      warning("%d observations were dropped because the corresponding probabilities for node %s were 0 or NaN.", dropped, NODE(i));

    if (debuglevel > 0)
      Rprintf("  > log-likelihood loss for node %s is %lf.\n", NODE(i), cur_loss);

    /* add the node contribution to the return value. */
    result += cur_loss;

  }/*FOR*/

  Free1D(columns);
  if (!c_is(fitted, "bn.fit.gnet"))
    Free1D(configs);

  UNPROTECT(3);
  return result;

}/*C_ENTROPY_LOSS*/

/* Gaussian loss for a single node. */
double c_gloss(int *cur, SEXP cur_parents, double *coefs, double *sd,
    void **columns, SEXP nodes, int ndata, double *per_sample,
    int allow_singular) {

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
      mean += ((double *)columns[p[j] - 1])[i] * coefs[j + 1];

    /* compute the log-likelihood of this observation. */
    if ((*sd < MACHINE_TOL) && !allow_singular)
      logprob = dnorm(((double *)columns[*cur])[i], mean, MACHINE_TOL, TRUE);
    else
      logprob = dnorm(((double *)columns[*cur])[i], mean, *sd, TRUE);

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
    SEXP data, SEXP nodes, int ndata, int nlevels, double *per_sample,
    int *dropped) {

int i = 0, *obs = NULL;
double logprob = 0, result = 0;
SEXP temp_df;

  /* get the target variable. */
  obs = INTEGER(VECTOR_ELT(data, *cur));
  /* get the parents' configurations. */
  if (length(cur_parents) > 0) {

    PROTECT(temp_df = c_dataframe_column(data, cur_parents, FALSE, FALSE));
    cfg(temp_df, configs, NULL);

    for (i = 0; i < ndata; i++) {

      /* propagate missing values from the parents configuration. */
      if (configs[i] == NA_INTEGER)
        logprob = NA_REAL;
      else
        logprob = log(prob[CMC(obs[i] - 1, configs[i], nlevels)]);

      UPDATE_LOSS(!R_FINITE(logprob) || ISNAN(logprob));

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else {

    for (i = 0; i < ndata; i++) {

      logprob = log(prob[obs[i] - 1]);

      UPDATE_LOSS(!R_FINITE(logprob) || ISNAN(logprob));

    }/*FOR*/

  }/*ELSE*/

  /* switch to the negentropy. */
  if (ndata > *dropped)
    result /= -(ndata - *dropped);
  else
    result = NA_REAL;

  return result;

}/*C_DLOSS*/

/* conditional Gaussian loss for a single node. */
double c_cgloss(int *cur, SEXP cur_parents, SEXP dparents, SEXP gparents,
    SEXP dlevels, double *coefs, double *sd, void **columns, SEXP nodes,
    int ndata, double *per_sample, int allow_singular, int *dropped) {

int i = 0, j = 0, *p = NULL, nparents = length(cur_parents);
int **dpar = NULL, *config = NULL;
int *ddp = INTEGER(dparents), *ggp = INTEGER(gparents);
int ndpar = length(dparents), ngpar = length(gparents), *nlvls = NULL;
double mean = 0, logprob = 0, result = 0, **gpar = NULL, *coefs_offset = NULL;
SEXP try;

  /* there is always at least one discrete parent. */
  PROTECT(try = match(nodes, cur_parents, 0));
  p = INTEGER(try);
  /* generate the discrete parents' configurations. */
  if (nparents == 1) {

    nparents = 0;
    config = (int *)columns[p[ddp[0] - 1] - 1];

  }/*THEN*/
  else {

    dpar = Calloc1D(ndpar, sizeof(int *));
    for (i = 0; i < ndpar; i++)
      dpar[i] = columns[p[ddp[i] - 1] - 1];
    nlvls = Calloc1D(ndpar, sizeof(int));
    for (i = 0; i < ndpar; i++)
      nlvls[i] = length(VECTOR_ELT(dlevels, i));
    config = Calloc1D(ndata, sizeof(int));
    c_fast_config(dpar, ndata, ndpar, nlvls, config, NULL, 1);

  }/*ELSE*/
  /* extract the continuous parents. */
  if (ngpar > 0) {

    gpar = Calloc1D(ngpar, sizeof(double *));
    for (i = 0; i < ngpar; i++)
      gpar[i] = columns[p[ggp[i] - 1] - 1];

  }/*THEN*/

  for (i = 0; i < ndata; i++) {

    if (config[i] == NA_INTEGER) {

      logprob = NA_REAL;

    }/*THEN*/
    else {

      coefs_offset = coefs + (ngpar + 1) * (config[i] - 1);

      /* compute the mean value for this observation. */
      mean = coefs_offset[0];

      for (j = 0; j < ngpar; j++)
        mean += gpar[j][i] * coefs_offset[j + 1];

      /* compute the log-likelihood of this observation. */
      if ((*(sd + config[i] - 1) < MACHINE_TOL) && !allow_singular)
        logprob = dnorm(((double *)columns[*cur])[i], mean, MACHINE_TOL, TRUE);
      else
        logprob = dnorm(((double *)columns[*cur])[i], mean, *(sd + config[i] - 1), TRUE);

    }/*ELSE*/

    UPDATE_LOSS((!R_FINITE(logprob) && !allow_singular) || ISNAN(logprob));

  }/*FOR*/

  UNPROTECT(1);

  if (ngpar > 0)
    Free1D(gpar);

  if (dpar) {

    Free1D(config);
    Free1D(nlvls);
    Free1D(dpar);

  }/*THEN*/

  /* switch to the negentropy. */
  if (ndata > *dropped)
    result /= -(ndata - *dropped);
  else
    result = NA_REAL;

  return result;

}/*C_CGLOSS*/

/* classification error of a single node as a loss function. */
SEXP class_err(SEXP reference, SEXP predicted) {

int i = 0, dropped = 0, ndata = length(reference);
int *r = INTEGER(reference), *p = INTEGER(predicted);
double err = 0;

  /* count how many elements differ (this assumes the levels of the two factors
   * are the same and in the same order); NAs are dropped. */
  for (i = 0; i < ndata; i++) {

    if ((r[i] == NA_INTEGER) || (p[i] == NA_INTEGER))
      dropped++;
    else if (r[i] != p[i])
      err++;

  }/*FOR*/

  /* rescale into a probability. */
  if (ndata > dropped)
    err /= (ndata - dropped);
  else
    err = NA_REAL;

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because of missing values.", dropped);

  return ScalarReal(err);

}/*CLASS_ERR*/
