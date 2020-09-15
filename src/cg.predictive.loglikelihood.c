#include "include/rcore.h"
#include "include/sets.h"
#include "include/scores.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/blas.h"

static double c_pcgnode(double *xx, double *xx2, double **gp, double **gp2,
    int ngp, int nobs, int nobs2, int *config, int *config2, int nconfig) {

int i = 0, j = 0;
double res = 0, fitted = 0, *beta = NULL, *sd = NULL;

  /* allocate the fitted values and the standard error. */
  sd = Calloc1D((!config) ? 1 : nconfig, sizeof(double));

  /* if the regression is conditional on config, iterate over its values,
   * otherwise fit using the whole sample. */
  if (!config) {

    beta = Calloc1D(ngp + 1, sizeof(double));

    c_ols(gp, xx, nobs, ngp, NULL, NULL, beta, sd, FALSE);

    /* compute the log-likelihood (singular models haze zero density). */
    if (*sd < MACHINE_TOL)
      res = R_NegInf;
    else
      for (j = 0; j < nobs2; j++) {

        for (i = 1, fitted = beta[0]; i < ngp + 1; i++)
          fitted += beta[j] * gp2[i - 1][j];

        res += dnorm(xx2[j], fitted, *sd, TRUE);

      }/*FOR*/

  }/*THEN*/
  else {

    beta = Calloc1D((ngp + 1) * nconfig, sizeof(double));

    c_cls(gp, xx, config, nobs, ngp, nconfig, NULL, NULL, beta, sd, FALSE);

    /* if the model is not singular compute the log-likelihood. */
    for (j = 0; j < nobs2; j++) {

      if (sd[config2[j] - 1] < MACHINE_TOL) {

        res = R_NegInf;
        goto end;

      }/*THEN*/

      for (i = 1, fitted = beta[0 + (config2[j] - 1) * (ngp + 1)]; i < ngp + 1; i++)
         fitted += beta[i + (config2[j] - 1) * (ngp + 1)] * gp2[i - 1][j];

      res += dnorm(xx2[j], fitted, sd[config2[j] - 1], TRUE);

    }/*FOR*/

  }/*ELSE*/

end:

  Free1D(beta);
  Free1D(sd);

  return res;

}/*C_PCGNODE*/

static double pcgnode(SEXP x, SEXP parents, SEXP x2, SEXP parents2, int *type,
    int np, int ndp, double *nparams) {

int i = 0, j = 0, k = 0, nobs = length(x), nobs2 = length(x2);
int nconfig = 0, ngp = np - ndp;
int **dp = NULL, **dp2 = NULL, *nlevels = NULL, *config = NULL, *config2 = NULL;
double res = 0, **gp = NULL, **gp2 = NULL;
SEXP temp;

  /* extract the discrete parents and build configurations. */
  dp = Calloc1D(ndp, sizeof(int *));
  dp2 = Calloc1D(ndp, sizeof(int *));
  gp = Calloc1D(ngp, sizeof(double *));
  gp2 = Calloc1D(ngp, sizeof(double *));
  config = Calloc1D(nobs, sizeof(int));
  config2 = Calloc1D(nobs, sizeof(int));
  nlevels = Calloc1D(ndp, sizeof(int));
  for (i = 0, j = 0, k = 0; i < np; i++)
    if (type[i] == INTSXP) {

      dp2[j] = INTEGER(VECTOR_ELT(parents2, i));
      temp = VECTOR_ELT(parents, i);
      dp[j] = INTEGER(temp);
      nlevels[j++] = NLEVELS(temp);

    }/*THEN*/
    else {

      gp2[k] = REAL(VECTOR_ELT(parents2, i));
      gp[k++] = REAL(VECTOR_ELT(parents, i));

    }/*ELSE*/

  c_fast_config(dp, nobs, ndp, nlevels, config, &nconfig, 1);
  c_fast_config(dp2, nobs2, ndp, nlevels, config2, NULL, 1);

  res = c_pcgnode(REAL(x), REAL(x2), gp, gp2, ngp, nobs, nobs2, config, config2,
          nconfig);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = nconfig * (ngp + 1);

  Free1D(dp);
  Free1D(dp2);
  Free1D(gp);
  Free1D(gp2);
  Free1D(config);
  Free1D(config2);
  Free1D(nlevels);

  return res;

}/*PCGNODE*/

double predictive_loglik_cgnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, int debugging) {

double loglik = 0;
int i = 0, nparents = 0, *type = NULL, cur_type = 0, dparents = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, new_t;
SEXP parent_vars, new_parents, config, config2;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  nparents = length(parents);
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
  PROTECT(new_t = c_dataframe_column(newdata, target, TRUE, FALSE));

  if (nparents == 0) {

    /* no parents, reuse the marginal likelihoods. */
    if (TYPEOF(data_t) == INTSXP)
      loglik = pdnode(data_t, new_t, nparams);
    else
      loglik = pgnode(data_t, new_t, nparams);

  }/*THEN*/
  else {

    /* extract the parents' columns from the data frame. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(new_parents = c_dataframe_column(newdata, parents, FALSE, FALSE));
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
        PROTECT(config2 = c_configurations(new_parents, TRUE, TRUE));
        loglik = cpdnode(data_t, config, new_t, config2, nparams);
        UNPROTECT(2);

      }/*ELSE*/

    }/*THEN*/
    else {

      if (dparents == 0)
        loglik = cpgnode(data_t, new_t, data, newdata, parents, nparams);
      else
        loglik = pcgnode(data_t, parent_vars, new_t, new_parents, type,
                   nparents, dparents, nparams);

    }/*ELSE*/

    Free1D(type);

    UNPROTECT(2);

  }/*ELSE*/

  if (debugging > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  UNPROTECT(2);

  return loglik;

}/*PREDICTIVE_LOGLIK_CGNODE*/
