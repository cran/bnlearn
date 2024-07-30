#include "../../include/rcore.h"
#include "../../fitted/fitted.h"
#include "../../core/data.table.h"
#include "../../core/allocations.h"
#include "loglikelihood.h"
#include "../../include/globals.h"


/* log-likelihood of individual observations for a Gaussian network. */
void bysample_gaussian_loglikelihood(fitted_bn bn, gdata dt, double *loglik,
    bool robust, bool debugging) {

int *pars = NULL;
double *obs = NULL, *coefs = NULL, sd = 0, *scratch = NULL;

  /* allocate the scratch space used to compute the expected value of the
   * distribution of each observation. */
  scratch = Calloc1D(dt.m.nobs, sizeof(double));

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... extract the relevant quantities from the data and the network... */
    obs = dt.col[i];
    coefs = bn.ldists[i].g.coefs;
    pars = bn.ldists[i].parents;
    sd = bn.ldists[i].g.sd;

    /* ... move away from singular distribtions ... */
    if ((sd < MACHINE_TOL) && robust)
      sd = MACHINE_TOL;

    /* ... initialize the log-likelihoods with the intercept... */
    for (int j = 0; j < dt.m.nobs; j++)
      scratch[j] = coefs[0];

    /* ... then add the effects of the parents one by one to compute the
     * expected values (no need to special case missing values, they propagate
     * correctly)... */
    for (int k = 0; k < bn.ldists[i].nparents; k++)
      for (int j = 0; j < dt.m.nobs; j++)
        scratch[j] += dt.col[pars[k]][j] * coefs[k + 1];

    /* ... and use them together with the standard error to compute the
     * log-likelihood. */
    for (int j = 0; j < dt.m.nobs; j++)
      loglik[j] += dnorm(obs[j], scratch[j], sd, TRUE);

  }/*FOR*/

  Free1D(scratch);

}/*BYSAMPLE_GAUSSIAN_LOGLIKELIHOOD*/

/* log-likelihood of a whole sample for a Gaussian network. */
double data_gaussian_loglikelihood(fitted_bn bn, gdata dt, double *scratch,
    bool propagate, bool loss, bool debugging) {

int ncomplete = 0, *pars = NULL;
double loglik = 0, node_loglik = 0;
double *obs = NULL, *coefs = NULL, sd = 0;
bool early_return = FALSE;

  /* if the data contain missing values for the nodes we are considering, and
   * we propagate them, the log-likelihood is necessarily NA. */
  if (propagate && check_locally_incomplete_data(bn, dt.m, debugging))
    return NA_REAL;

  /* if any of the coefficients of the nodes we are considering is NA, or the
   * standard error is NA, then the log-likelihood is NA. */
  for (int i = 0; i < bn.nnodes; i++) {

    if (!dt.m.flag[i].fixed)
      continue;

    if (ISNAN(bn.ldists[i].g.sd)) {

      early_return = TRUE;
      goto unidentifiable_model;

    }/*THEN*/

    for (int j = 0; j < bn.ldists[i].g.ncoefs; j++)
      if (ISNAN(bn.ldists[i].g.coefs[j])) {

        early_return = TRUE;
        goto unidentifiable_model;

      }/*THEN*/

unidentifiable_model:
    if (early_return) {

      if (debugging)
        Rprintf("* unidentifiable model in node %s, the log-likelihood is NA.\n",
            bn.labels[i]);
      return NA_REAL;

    }/*THEN*/

  }/*FOR*/

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging && !loss)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... reset the log-likelihood accumulator... */
    node_loglik = 0;
    ncomplete = 0;

    /* ... extract the relevant quantities from the data and the network...*/
    obs = dt.col[i];
    coefs = bn.ldists[i].g.coefs;
    pars = bn.ldists[i].parents;
    sd = bn.ldists[i].g.sd;

    /* ... initialize the log-likelihoods with the intercept... */
    for (int j = 0; j < dt.m.nobs; j++)
      scratch[j] = coefs[0];

    /* ... then add the effects of the parents one by one to compute the
     * expected values (no need to special case missing values, they propagate
     * correctly)... */
    for (int k = 0; k < bn.ldists[i].nparents; k++)
      for (int j = 0; j < dt.m.nobs; j++)
        scratch[j] += dt.col[pars[k]][j] * coefs[k + 1];

    /* ... and use them together with the standard error to compute the
     * log-likelihood. */
    for (int j = 0; j < dt.m.nobs; j++) {

      if (ISNAN(obs[j]) || ISNAN(scratch[j]))
        continue;

      node_loglik += dnorm(obs[j], scratch[j], sd, TRUE);
      ncomplete++;

    }/*FOR*/

    /* scale the likelihood to compensate for any missing values (which will
     * not be propagated as a result), or return -Inf if there are no
     * locally-complete observations. */
    if (ncomplete == 0)
      node_loglik = R_NegInf;
    else if (ncomplete < dt.m.nobs)
      node_loglik = node_loglik / ncomplete * dt.m.nobs;

    if (loss) {

      if (debugging)
        Rprintf("  > log-likelihood loss for node %s is %lf.\n",
          bn.labels[i], - node_loglik / dt.m.nobs);

    }/*THEN*/
    else {

      if (debugging) {

        Rprintf("  > %d locally-complete observations out of %d.\n",
          ncomplete, dt.m.nobs);
        Rprintf("  > log-likelihood is %lf.\n", node_loglik);

      }/*THEN*/

    }/*ELSE*/

    /* cumulate the log-likelihood. */
    loglik += node_loglik;
    /* if the log-likelihood is NA or -Inf it will never change value again. */
    if (ISNAN(loglik) || (loglik == R_NegInf))
      break;

  }/*FOR*/

  return loglik;

}/*DATA_GAUSSIAN_LOGLIKELIHOOD*/

