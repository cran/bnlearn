#include "include/rcore.h"
#include "include/sets.h"
#include "include/scores.h"
#include "include/dataframe.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/blas.h"

double ccgloglik(SEXP x, SEXP parents, int *type, int np, int ndp,
    double *nparams) {

int i = 0, j = 0, k = 0, nobs = length(x), nconfig = 0, ngp = np - ndp;
int **dp = NULL, *nlevels = NULL, *config = NULL;
double res = 0, **gp = NULL;
SEXP temp;

  /* extract the discrete parents and build configurations. */
  dp = Calloc1D(ndp, sizeof(int *));
  gp = Calloc1D(ngp, sizeof(double *));
  config = Calloc1D(nobs, sizeof(int));
  nlevels = Calloc1D(ndp, sizeof(int));
  for (i = 0, j = 0, k = 0; i < np; i++)
    if (type[i] == INTSXP) {

      temp = VECTOR_ELT(parents, i);
      dp[j] = INTEGER(temp);
      nlevels[j++] = NLEVELS(temp);

    }/*THEN*/
    else {

      gp[k++] = REAL(VECTOR_ELT(parents, i));

    }/*ELSE*/
  c_fast_config(dp, nobs, ndp, nlevels, config, &nconfig, 1);

  res = c_fast_ccgloglik(REAL(x), gp, ngp, nobs, config, nconfig);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = nconfig * (ngp + 1);

  Free1D(dp);
  Free1D(gp);
  Free1D(config);
  Free1D(nlevels);

  return res;

}/*CCGLOGLIK*/

double c_fast_ccgloglik(double *xx, double **gp, int ngp, int nobs, int *config,
    int nconfig) {

int i = 0, j = 0, k = 0, l = 0, batch = 0;
double res = 0, *subset = NULL, *fitted = NULL, *response = NULL;
long double sd = 0;

  /* if the regression is conditional on config, iterate over its values,
   * otherwise fit using the whole sample. */
  if (!config) {

    /* allocate the array of continuous columns (plus the intercept). */
    subset = Calloc1D(nobs * (ngp + 1), sizeof(double));
    /* allocate the fitted values and the response variable. */
    fitted = Calloc1D(nobs, sizeof(double));
    /* initialize the intercept. */
    for (j = 0; j < nobs; j++)
      subset[CMC(j, 0, batch)] = 1;
    /* copy the data from the continuous parents. */
    for (k = 0; k < ngp; k++)
      memcpy(subset + (k + 1) * nobs, gp[k], nobs * sizeof(double));

    /* estimate OLS using the QR decomposition. */
    c_qr_ols(subset, xx, nobs, ngp + 1, fitted, &sd);
    /* compute the log-likelihood (singular models haze zero density). */
    if (sd < MACHINE_TOL)
      res = R_NegInf;
    else
      for (j = 0; j < nobs; j++)
        res += dnorm(xx[j], fitted[j], (double)sd, TRUE);

    Free1D(subset);
    Free1D(fitted);

  }/*THEN*/
  else {

    /* iterate over the configurations. */
    for (i = 1; i <= nconfig; i++) {

      /* count how many observations have this parents' configuration. */
      for (j = 0, batch = 0; j < nobs; j++)
        if (config[j] == i)
          batch++;

      /* no observations for this configuration, skip. */
      if (batch == 0)
        continue;

      /* allocate the array of continuous columns (plus the intercept). */
      subset = Calloc1D(batch * (ngp + 1), sizeof(double));
      /* allocate the fitted values and the response variable. */
      fitted = Calloc1D(batch, sizeof(double));
      response = Calloc1D(batch, sizeof(double));
      /* initialize the intercept. */
      for (j = 0; j < batch; j++)
        subset[CMC(j, 0, batch)] = 1;
      /* copy the data from the continuous parents. */
      for (k = 0; k < ngp; k++)
        for (j = 0, l = 0; j < nobs; j++)
          if (config[j] == i)
            subset[CMC(l++, k + 1, batch)] = gp[k][j];
      /* do the same with the response variable. */
      for (j = 0, l = 0; j < nobs; j++)
        if (config[j] == i)
           response[l++] = xx[j];

      /* estimate OLS using the QR decomposition. */
      c_qr_ols(subset, response, batch, ngp + 1, fitted, &sd);
      /* compute the log-likelihood (singular models haze zero density). */
      if (sd < MACHINE_TOL)
        res = R_NegInf;
      else
        for (j = 0; j < batch; j++)
          res += dnorm(response[j], fitted[j], (double)sd, TRUE);

      Free1D(subset);
      Free1D(fitted);
      Free1D(response);

      /* if the log-likelihood has degenerated into a point mass, there's no
       * reason to go on with more batches. */
      if (!R_FINITE(res))
        break;

    }/*FOR*/

  }/*ELSE*/

  return res;

}/*CCGLOGLIK*/

double loglik_cgnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel) {

double loglik = 0;
int i = 0, nparents = 0, *type = NULL, cur_type = 0, dparents = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, parent_vars, config;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  nparents = length(parents);
  /* extract the node's column from the data frame. */
  data_t = c_dataframe_column(data, target, TRUE, FALSE);

  if (nparents == 0) {

    /* no parents, reuse the marginal likelihoods. */
    if (TYPEOF(data_t) == INTSXP)
      loglik = dlik(data_t, nparams);
    else
      loglik = glik(data_t, nparams);

  }/*THEN*/
  else {

    /* extract the parents' columns from the data frame. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    /* find out which kind of parents the node has, and count discrete parents
     * (all the others are continuous). */
    type = Calloc1D(nparents, sizeof(int));

    for (i = 0; i < nparents; i++) {

      cur_type = TYPEOF(VECTOR_ELT(parent_vars, i));
      dparents += (cur_type == INTSXP);
      type[i] = cur_type;

    }/*FOR*/

    if (TYPEOF(data_t) == INTSXP) {

      /* discrete nodes are not allowed to have continuous parents; then reuse
       * the conditional discrete likelihood.*/
     if (dparents != nparents) {

       loglik = R_NegInf;

     }/*THEN*/
     else {

       PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
       loglik = cdlik(data_t, config, nparams);
       UNPROTECT(1);

     }/*ELSE*/

    }/*THEN*/
    else {

      if (dparents == 0)
        loglik = cglik(data_t, data, parents, nparams);
      else
        loglik = ccgloglik(data_t, parent_vars, type, nparents, dparents,
                   nparams);

    }/*ELSE*/

    Free1D(type);

    UNPROTECT(1);

  }/*ELSE*/

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_CGNODE*/

