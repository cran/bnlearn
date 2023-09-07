#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/sets.h"
#include "scores.h"
#include "../minimal/data.frame.h"
#include "../minimal/common.h"
#include "../include/globals.h"
#include "../math/linear.algebra.h"

double c_fast_ccgnal(double *xx, double **gp, int ngp, int nobs,
    int *config, int nconfig, double k) {

int i = 0, j = 0, nc = 0;
double res = 0, *fitted = NULL, *sd = NULL;

  /* allocate the fitted values and the standard error. */
  fitted = Calloc1D(nobs, sizeof(double));
  sd = Calloc1D((!config) ? 1 : nconfig, sizeof(double));

  /* if the regression is conditional on config, iterate over its values,
   * otherwise fit using the whole sample. */
  if (!config) {

    c_ols(gp, xx, nobs, ngp, fitted, NULL, NULL, sd, &nc, TRUE);

    /* compute the log-likelihood. */
    if ((*sd < MACHINE_TOL) || (nc == 0)) {

      /* singular distributions and those with no locally-complete observations
       * are scored to prevent their selection in structure learning. */
      res = R_NegInf;

    }/*THEN*/
    else {

      for (j = 0; j < nobs; j++)
        if (!ISNAN(xx[j]) && !ISNAN(fitted[j]) && !ISNAN(*sd))
          res += dnorm(xx[j], fitted[j], *sd, TRUE);

    }/*ELSE*/

  }/*THEN*/
  else {

    c_cls(gp, xx, config, nobs, ngp, nconfig, fitted, NULL, NULL, sd, &nc,
      TRUE);

    /* if any standard error is zero, the model is singular; if it is undefined,
     * there are no locally-complete observations. */
    for (i = 0; i < nconfig; i++) {

      if ((sd[i] < MACHINE_TOL) || ISNAN(sd[i])) {

        res = R_NegInf;
        goto end;

      }/*FOR*/

    }/*FOR*/

    /* if the model is not singular compute the log-likelihood. */
    for (j = 0; j < nobs; j++)
      if (!ISNAN(xx[j]) && !ISNAN(fitted[j]) && !ISNAN(sd[config[j] - 1]))
        res += dnorm(xx[j], fitted[j], sd[config[j] - 1], TRUE);

  }/*ELSE*/

  /* average over the locally-complete data... */
  res /= nc;
  /* and add the penalty term, scaled by the original sample size. */
  res -= k / nobs * (ngp + 2) * ((!config) ? 1 : nconfig);

end:

  Free1D(fitted);
  Free1D(sd);

  return res;

}/*C_FAST_CCGNAL*/

double ccgnal(SEXP x, SEXP parents, int *type, int np, int ndp, double k) {

int i = 0, j = 0, l = 0, nobs = length(x), nconfig = 0, ngp = np - ndp;
int **dp = NULL, *nlevels = NULL, *config = NULL;
double res = 0, **gp = NULL;
SEXP temp;

  /* extract the discrete parents and build configurations. */
  dp = Calloc1D(ndp, sizeof(int *));
  gp = Calloc1D(ngp, sizeof(double *));
  config = Calloc1D(nobs, sizeof(int));
  nlevels = Calloc1D(ndp, sizeof(int));
  for (i = 0, j = 0, l = 0; i < np; i++)
    if (type[i] == INTSXP) {

      temp = VECTOR_ELT(parents, i);
      dp[j] = INTEGER(temp);
      nlevels[j++] = NLEVELS(temp);

    }/*THEN*/
    else {

      gp[l++] = REAL(VECTOR_ELT(parents, i));

    }/*ELSE*/
  c_fast_config(dp, nobs, ndp, nlevels, config, &nconfig, 1);

  res = c_fast_ccgnal(REAL(x), gp, ngp, nobs, config, nconfig, k);

  Free1D(dp);
  Free1D(gp);
  Free1D(config);
  Free1D(nlevels);

  return res;

}/*CCGNAL*/

double nal_cgnode(SEXP target, SEXP x, SEXP data, double k, bool debugging) {

double nal = 0;
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
      nal = nal_dnode_root(data_t, k);
    else
      nal = glik_incomplete(data_t, k);

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

       nal = R_NegInf;

     }/*THEN*/
     else {

       PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
       nal = nal_dnode_parents(data_t, config, k);
       UNPROTECT(1);

     }/*ELSE*/

    }/*THEN*/
    else {

      if (dparents == 0)
        nal = cglik_incomplete(data_t, data, parents, k);
      else
        nal = ccgnal(data_t, parent_vars, type, nparents, dparents, k);

    }/*ELSE*/

    Free1D(type);

    UNPROTECT(1);

  }/*ELSE*/

  if (debugging)
    Rprintf("  > log-likelihood is %lf.\n", nal);

  UNPROTECT(1);

  return nal;

}/*NAL_CGNODE*/

