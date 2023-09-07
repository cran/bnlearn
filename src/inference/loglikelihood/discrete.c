#include "../../include/rcore.h"
#include "../../fitted/fitted.h"
#include "../../core/data.table.h"
#include "../../core/allocations.h"
#include "../../core/contingency.tables.h"
#include "../../core/math.functions.h"
#include "../../core/sets.h"
#include "../../math/linear.algebra.h"
#include "loglikelihood.h"

/* log-likelihood of individual observations for a discrete network. */
void bysample_discrete_loglikelihood(fitted_bn bn, ddata dt, double *loglik,
    bool debugging) {

int *obs = NULL, *parcfgs = NULL;
double *cpt = NULL;
bool locally_complete = FALSE;
ddata local_data = { 0 };

  /* allocate a second data table to hold the data that are local to the node's
   * parents. */
  local_data = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* parent configurations to index the parameters in the probability tables. */
  parcfgs = Calloc1D(dt.m.nobs, sizeof(int));

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... get its values from the data and the (conditional) probability table
     * from the network... */
    obs = dt.col[i];
    cpt = bn.ldists[i].d.cpt;
    locally_complete = dt.m.flag[i].complete;

    if (bn.ldists[i].nparents == 0) {

      /* ... if the node is a root node, the value of the observation is the
       * index in the probability table ... */
      if (locally_complete) {

        for (int j = 0; j < dt.m.nobs; j++)
          loglik[j] += log(cpt[obs[j] - 1]);

      }/*THEN*/
      else {

        for (int j = 0; j < dt.m.nobs; j++) {

          if (obs[j] == NA_INTEGER)
            loglik[j] = NA_REAL;
           else
            loglik[j] += log(cpt[obs[j] - 1]);

        }/*FOR*/

      }/*ELSE*/

    }/*THEN*/
    else {

      /* ... if the node has parents, extract them from the data... */
      ddata_subset_columns(&dt, &local_data, bn.ldists[i].parents,
         bn.ldists[i].nparents);
      /* ... compute their configurations... */
      c_fast_config(local_data.col, local_data.m.nobs, local_data.m.ncols,
          local_data.nlvl, parcfgs, NULL, 0);

      for (int k = 0; k < local_data.m.ncols; k++)
        locally_complete &= local_data.m.flag[k].complete;

      /* ... and the two together index the conditional probability table. */
      if (locally_complete) {

        for (int j = 0; j < dt.m.nobs; j++)
          loglik[j] += log(cpt[CMC(obs[j] - 1, parcfgs[j],
                           bn.ldists[i].d.dims[0])]);

      }/*THEN*/
      else {

        for (int j = 0; j < dt.m.nobs; j++) {

          if ((obs[j] == NA_INTEGER) || (parcfgs[j] == NA_INTEGER))
            loglik[j] = NA_REAL;
          else
            loglik[j] += log(cpt[CMC(obs[j] - 1, parcfgs[j],
                             bn.ldists[i].d.dims[0])]);

        }/*FOR*/

      }/*ELSE*/

    }/*ELSE*/

  }/*FOR*/

  Free1D(parcfgs);
  FreeDDT(local_data);

}/*BYSAMPLE_DISCRETE_LOGLIKELIHOOD*/

/* log-likelihood of a whole sample for a discrete network. */
double data_discrete_loglikelihood(fitted_bn bn, ddata dt, bool propagate,
    bool debugging) {

int max_nlvls = 0, cumdim = 0, max_cfgs = 0, *parcfgs = NULL;
double loglik = 0, node_loglik = 0;
counts1d freq = { 0 };
counts2d freq2 = { 0 };
ddata local_data = { 0 };

  /* if the data contain missing values for the nodes we are considering, and
   * we propagate them, the log-likelihood is necessarily NA. */
  if (propagate && check_locally_incomplete_data(bn, dt.m, debugging))
    return NA_REAL;

  /* allocate a second data table to hold the data that are local to the node's
   * parents. */
  local_data = empty_ddata(dt.m.nobs, dt.m.ncols);
  /* parent configurations to index the parameters in the probability tables. */
  parcfgs = Calloc1D(dt.m.nobs, sizeof(int));
  /* find out the largest number of levels that a variable can take to bound the
   * size of the contingency tables for the counts. */
  max_nlvls = dt.nlvl[i_which_max(dt.nlvl, dt.m.ncols) - 1];
  for (int i = 0; i < bn.nnodes; i++) {

    cumdim = 1;
    for (int j = 1; j < bn.ldists[i].d.ndims; j++)
      cumdim *= bn.ldists[i].d.dims[j];

    max_cfgs = (max_cfgs > cumdim) ? max_cfgs : cumdim;

  }/*FOR*/

  freq = new_1d_table(max_nlvls);
  freq2 = new_2d_table(max_nlvls, max_cfgs, FALSE);

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... reset the log-likelihood accumulator... */
    node_loglik = 0;

    if (bn.ldists[i].nparents == 0) {

      /* ... compute the frequencies... */
      resize_1d_table(dt.nlvl[i], &freq);
      refill_1d_table(dt.col[i], &freq, dt.m.nobs);
      /* ... if there are frequencies... */
      if (freq.nobs == 0)
        node_loglik = R_NegInf;
      else {

        /* ... combine them with the probabilities to compute the
         * log-likelihood... */
        for (int j = 0; j < freq.llx; j++)
          if (freq.n[j] > 0)
            node_loglik += freq.n[j] * log(bn.ldists[i].d.cpt[j]);
        /* ... and scale the log-likelihood to compensate for any missing values
         * (which will not be propagated as a result). */
        if (freq.nobs < dt.m.nobs)
          node_loglik = node_loglik / freq.nobs * dt.m.nobs;

      }/*ELSE*/

    }/*THEN*/
    else {

      /* ... construct the parent configurations (with a 1-offset to make
       * fill_2d_table() happy)... */
      ddata_subset_columns(&dt, &local_data, bn.ldists[i].parents,
         bn.ldists[i].nparents);
      c_fast_config(local_data.col, local_data.m.nobs, local_data.m.ncols,
          local_data.nlvl, parcfgs, &cumdim, 1);
      /* ... compute the frequencies... */
      resize_2d_table(dt.nlvl[i], cumdim, &freq2);
      refill_2d_table(dt.col[i], parcfgs, &freq2, dt.m.nobs);
      /* ... if there are usable locally-complete observations... */
      if (freq2.nobs == 0)
        node_loglik = R_NegInf;
      else {

        /* ... combine them with the probabilities to compute the
         * log-likelihood... */
        for (int j = 0; j < freq2.llx; j++)
          for (int k = 0; k < freq2.lly; k++)
            if (freq2.n[j][k] > 0)
              node_loglik += freq2.n[j][k] *
                       log(bn.ldists[i].d.cpt[CMC(j, k, bn.ldists[i].d.dims[0])]);
        /* ... and scale the log-likelihood to compensate for any missing values
         * (which will not be propagated as a result). */
        if (freq2.nobs < dt.m.nobs)
          node_loglik = node_loglik / freq2.nobs * dt.m.nobs;

      }/*ELSE*/

    }/*ELSE*/

    if (debugging) {

      Rprintf("  > %d locally-complete observations out of %d.\n",
        (bn.ldists[i].nparents == 0) ? freq.nobs: freq2.nobs, dt.m.nobs);
      Rprintf("  > log-likelihood is %lf.\n", node_loglik);

    }/*THEN*/

    /* cumulate the log-likelihood. */
    loglik += node_loglik;
    /* if the log-likelihood is NA or -Inf it will never change value again. */
    if (ISNAN(loglik) || (loglik == R_NegInf))
      break;

  }/*FOR*/

  Free1D(parcfgs);
  resize_1d_table(max_nlvls, &freq);
  Free1DTAB(freq);
  resize_2d_table(max_nlvls, max_cfgs, &freq2);
  Free2DTAB(freq2);
  FreeDDT(local_data);

  return loglik;

}/*DATA_DISCRETE_LOGLIKELIHOOD*/

