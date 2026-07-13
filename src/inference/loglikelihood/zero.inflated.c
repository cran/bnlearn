#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/data.table.h"
#include "../../fitted/fitted.h"
#include "../../math/hyperpoisson.h"
#include "../../math/xnegbin.h"
#include "../../scores/scores.h"
#include "loglikelihood.h"

/* log-likelihood of individual observations for a zero-inflated network. */
void bysample_zeroinflated_loglikelihood(fitted_bn bn, tabular dt, double *loglik,
    bool robust, bool debugging) {

int *pars = NULL, npars = 0;
double *obs = NULL, *scratch1 = NULL, *scratch2 = NULL;
tabular sub = { 0 };

  /* allocate a second data set for storing the parents. */
  sub = empty_tabular(dt.m.nobs, 0, dt.m.ncols);
  /* allocate the scratch space for zero-inflation probability and parameter
   * values for each observation. */
  scratch1 = Calloc1D(dt.m.nobs, sizeof(double));
  scratch2 = Calloc1D(dt.m.nobs, sizeof(double));

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... extract the relevant quantities from the data and the network... */
    obs = dt.ccol[i];
    pars = bn.ldists[i].parents;
    npars = bn.ldists[i].nparents;

    /* prepare the current subset. */
    tabular_subset_columns(&dt, &sub, pars, npars);

    if (bn.node_types[i] == ZIHPNODE) {

      double *zinf_prob = scratch1, *intensity = scratch2;
      double dispersion = bn.ldists[i].zihp.dispersion;

      /* translate the regression coefficients into canonical parameters... */
      zihp_coefs_to_params(zinf_prob, bn.ldists[i].zihp.inflation,
          intensity, bn.ldists[i].zihp.intensity, sub);

      /* corner case: at dispersion exactly zero the rising factorial and the
       * confluent hypergeometric function are individually singular, but the
       * hyper-Poisson has a well-defined limit -- a shifted Poisson,
       * X - 1 ~ Poisson(intensity), with no mass at zero (which then comes only
       * from the zero-inflation component); this matches the sampler in
       * rhypois(). For any dispersion > 0, dhypois() below is finite and
       * already tends to this same limit (via the log-gamma rising factorial),
       * so no special case is needed there. */
      if (dispersion == 0) {

        for (int j = 0; j < dt.m.nobs; j++)
          if (obs[j] == 0)
            loglik[j] += log(zinf_prob[j]);
          else
            loglik[j] += log(1 - zinf_prob[j]) +
                           dpois(obs[j] - 1, intensity[j], TRUE);

        continue;

      }/*THEN*/

      /* ... and use them together with the dispersion to compute the
       * log-likelihood. */
      for (int j = 0; j < dt.m.nobs; j++) {

        if (obs[j] == 0)
          loglik[j] += log(zinf_prob[j] + (1 - zinf_prob[j]) *
                         dhypois(obs[j], intensity[j], dispersion, FALSE));
        else
          loglik[j] += log(1 - zinf_prob[j]) +
                         dhypois(obs[j], intensity[j], dispersion, TRUE);

      }/*FOR*/

    }/*THEN*/
    else if (bn.node_types[i] == ZINBNODE) {

      double *zinf_prob = scratch1, *succ_prob = scratch2, fail_num = 0;
      fail_num = bn.ldists[i].zinb.failures;

      /* translate the regression coefficients into canonical parameters. */
      zinb_coefs_to_params(zinf_prob, bn.ldists[i].zinb.inflation,
          succ_prob, bn.ldists[i].zinb.prsucc, sub);

      /* ... and use them together with the number of failures to compute the
       * log-likelihood. */
      for (int j = 0; j < dt.m.nobs; j++) {

        if (obs[j] == 0)
          loglik[j] += log(zinf_prob[j] + (1 - zinf_prob[j]) *
                         dxnegbin(obs[j], succ_prob[j], fail_num, FALSE));
        else
          loglik[j] += log(1 - zinf_prob[j]) +
                         dxnegbin(obs[j], succ_prob[j], fail_num, TRUE);

      }/*FOR*/

    }/*THEN*/

  }/*FOR*/

  FreeTAB(sub);
  Free1D(scratch1);
  Free1D(scratch2);

}/*BYSAMPLE_ZEROINFLATED_LOGLIKELIHOOD*/

/* log-likelihood of a whole sample for a zero-inflated network. */
double data_zeroinflated_loglikelihood(fitted_bn bn, tabular dt, double *scratch,
    bool propagate, bool debugging) {

int ncomplete = 0, *pars = NULL, npars = 0;
double loglik = 0, node_loglik = 0;
double *obs = NULL, *coefs = NULL, *scratch1 = NULL, *scratch2 = NULL;
bool early_return = FALSE;
tabular sub = { 0 };

  /* if the data contain missing values for the nodes we are considering, and
   * we propagate them, the log-likelihood is necessarily NA. */
  if (propagate && check_locally_incomplete_data(bn, dt.m, debugging))
    return NA_REAL;

  /* if any of the coefficients of the nodes we are considering is NA, or the
   * standard error is NA, then the log-likelihood is NA. */
  for (int i = 0; i < bn.nnodes; i++) {

    if (!dt.m.flag[i].fixed)
      continue;

    if (bn.node_types[i] == ZIHPNODE) {

      if (ISNAN(bn.ldists[i].zihp.dispersion)) {

        early_return = TRUE;
        goto unidentifiable_model;

      }/*THEN*/

      for (int j = 0; j < bn.ldists[i].zihp.ncoefs; j++)
        if (ISNAN(bn.ldists[i].zihp.inflation[j]) ||
            ISNAN(bn.ldists[i].zihp.intensity[j])) {

          early_return = TRUE;
          goto unidentifiable_model;

        }/*THEN*/
    }/*THEN*/
    else if (bn.node_types[i] == ZINBNODE) {

      if (ISNAN(bn.ldists[i].zinb.failures)) {

        early_return = TRUE;
        goto unidentifiable_model;

      }/*THEN*/

      for (int j = 0; j < bn.ldists[i].zinb.ncoefs; j++)
        if (ISNAN(bn.ldists[i].zinb.inflation[j]) ||
            ISNAN(bn.ldists[i].zinb.prsucc[j])) {

          early_return = TRUE;
          goto unidentifiable_model;

        }/*THEN*/

    }/*THEN*/

unidentifiable_model:
    if (early_return) {

      if (debugging)
        Rprintf("* unidentifiable model in node %s, the log-likelihood is NA.\n",
            bn.labels[i]);
      return NA_REAL;

    }/*THEN*/

  }/*FOR*/

  /* allocate a second data set for storing the parents. */
  sub = empty_tabular(dt.m.nobs, 0, dt.m.ncols);
  /* allocate the scratch space for zero-inflation probability and parameter
   * values for each observation. */
  scratch1 = Calloc1D(dt.m.nobs, sizeof(double));
  scratch2 = Calloc1D(dt.m.nobs, sizeof(double));

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... reset the log-likelihood accumulator... */
    node_loglik = 0;
    ncomplete = 0;

    /* ... extract the relevant quantities from the data and the network... */
    obs = dt.ccol[i];
    pars = bn.ldists[i].parents;
    npars = bn.ldists[i].nparents;

    /* prepare the current subset. */
    tabular_subset_columns(&dt, &sub, pars, npars);

    if (bn.node_types[i] == ZIHPNODE) {

      double *zinf_prob = scratch1, *intensity = scratch2;
      double dispersion = bn.ldists[i].zihp.dispersion;

      /* translate the regression coefficients into canonical parameters... */
      zihp_coefs_to_params(zinf_prob, bn.ldists[i].zihp.inflation,
          intensity, bn.ldists[i].zihp.intensity, sub);

      /* corner case: shifted Poisson at dispersion exactly zero, consistent
       * with the by-sample log-likelihood and the sampler in rhypois(); for any
       * dispersion > 0 dhypois() below is finite and tends to this limit. */
      if (dispersion == 0) {

        for (int j = 0; j < dt.m.nobs; j++) {

          if (ISNAN(obs[j]) || ISNAN(zinf_prob[j]) || ISNAN(intensity[j]))
            continue;

          if (obs[j] == 0)
            node_loglik += log(zinf_prob[j]);
          else
            node_loglik += log(1 - zinf_prob[j]) +
                             dpois(obs[j] - 1, intensity[j], TRUE);

          ncomplete++;

        }/*FOR*/

      }/*THEN*/
      else {

        /* ... and use them together with the dispersion to compute the
         * log-likelihood. */
        for (int j = 0; j < dt.m.nobs; j++) {

          if (ISNAN(obs[j]) || ISNAN(zinf_prob[j]) || ISNAN(intensity[j]) ||
              ISNAN(dispersion))
            continue;

          if (obs[j] == 0)
            node_loglik += log(zinf_prob[j] + (1 - zinf_prob[j]) *
                             dhypois(obs[j], intensity[j], dispersion, FALSE));
          else
            node_loglik += log(1 - zinf_prob[j]) +
                             dhypois(obs[j], intensity[j], dispersion, TRUE);

          ncomplete++;

        }/*FOR*/

      }/*ELSE*/

    }/*THEN*/
    else if (bn.node_types[i] == ZINBNODE) {

      double *zinf_prob = scratch1, *succ_prob = scratch2;
      double fail_num = bn.ldists[i].zinb.failures;

      /* translate the regression coefficients into canonical parameters. */
      zinb_coefs_to_params(zinf_prob, bn.ldists[i].zinb.inflation,
          succ_prob, bn.ldists[i].zinb.prsucc, sub);

      /* ... and use them together with the number of failures to compute the
       * log-likelihood. */
      for (int j = 0; j < dt.m.nobs; j++) {

        if (ISNAN(obs[j]) || ISNAN(zinf_prob[j]) || ISNAN(succ_prob[j]) ||
            ISNAN(fail_num))
          continue;

        if (obs[j] == 0)
          node_loglik += log(zinf_prob[j] + (1 - zinf_prob[j]) *
                           dxnegbin(obs[j], succ_prob[j], fail_num, FALSE));
        else
          node_loglik += log(1 - zinf_prob[j]) +
                           dxnegbin(obs[j], succ_prob[j], fail_num, TRUE);

        ncomplete++;

      }/*FOR*/

    }/*THEN*/

    /* scale the likelihood to compensate for any missing values (which will
     * not be propagated as a result), or return -Inf if there are no
     * locally-complete observations. */
    if (ncomplete == 0)
      node_loglik = NA_REAL;
    else if (ncomplete < dt.m.nobs)
      node_loglik = node_loglik / ncomplete * dt.m.nobs;

    if (debugging) {

      Rprintf("  > %d locally-complete observations out of %d.\n",
        ncomplete, dt.m.nobs);
      Rprintf("  > log-likelihood is %lf.\n", node_loglik);

    }/*THEN*/

    /* cumulate the log-likelihood. */
    loglik += node_loglik;
    /* if the log-likelihood is NA or -Inf it will never change value again. */
    if (ISNAN(loglik) || (loglik == R_NegInf))
      break;

  }/*FOR*/

  FreeTAB(sub);
  Free1D(scratch1);
  Free1D(scratch2);

  return loglik;

}/*DATA_ZEROINFLATED_LOGLIKELIHOOD*/

