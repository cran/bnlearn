#include "include/rcore.h"
#include "include/sets.h"
#include "include/scores.h"
#include "include/data.frame.h"
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

  /* we may want to store the number of parameters (one intercept, one
   * regression coefficient for each continuous parent, one residuals
   * standard error for each configuration of the discrete parents). */
  if (nparams)
    *nparams = nconfig * (ngp + 2);

  Free1D(dp);
  Free1D(gp);
  Free1D(config);
  Free1D(nlevels);

  return res;

}/*CCGLOGLIK*/

double c_fast_ccgloglik(double *xx, double **gp, int ngp, int nobs, int *config,
    int nconfig) {

int i = 0, j = 0;
double res = 0, *fitted = NULL, *sd = NULL;

  /* allocate the fitted values and the standard error. */
  fitted = Calloc1D(nobs, sizeof(double));
  sd = Calloc1D((!config) ? 1 : nconfig, sizeof(double));

  /* if the regression is conditional on config, iterate over its values,
   * otherwise fit using the whole sample. */
  if (!config) {

    c_ols(gp, xx, nobs, ngp, fitted, NULL, NULL, sd, FALSE);

    /* compute the log-likelihood (singular models haze zero density). */
    if (*sd < MACHINE_TOL)
      res = R_NegInf;
    else
      for (j = 0; j < nobs; j++)
        res += dnorm(xx[j], fitted[j], *sd, TRUE);

  }/*THEN*/
  else {

    c_cls(gp, xx, config, nobs, ngp, nconfig, fitted, NULL, NULL, sd, FALSE);

    /* if any standard error is zero, the model is singular and has density
     * zero. */
    for (i = 0; i < nconfig; i++) {

      if (sd[i] < MACHINE_TOL) {

        res = R_NegInf;
        goto end;

     }/*FOR*/

    }/*FOR*/

    /* if the model is not singular compute the log-likelihood. */
    for (j = 0; j < nobs; j++)
      res += dnorm(xx[j], fitted[j], sd[config[j] - 1], TRUE);

  }/*ELSE*/

end:

  Free1D(fitted);
  Free1D(sd);

  return res;

}/*CCGLOGLIK*/

double loglik_cgnode(SEXP target, SEXP x, SEXP data, double *nparams,
    bool debugging) {

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
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));

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

  if (debugging)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  UNPROTECT(1);

  return loglik;

}/*LOGLIK_CGNODE*/

