#include "../../include/rcore.h"
#include "../../fitted/fitted.h"
#include "../../core/data.table.h"
#include "../../core/allocations.h"
#include "../../core/contingency.tables.h"
#include "../../core/math.functions.h"
#include "../../core/sets.h"
#include "../../math/linear.algebra.h"
#include "loglikelihood.h"

/* log-likelihood of individual observations for a conditional Gaussian
 * network. */
void bysample_clgaussian_loglikelihood(fitted_bn bn, cgdata dt, double *loglik,
    bool debugging) {

int *pars = NULL, *dobs = NULL, *parcfgs = NULL, ncoefs = 0;
double *gobs = NULL, *coefs = NULL, sd = 0, *sds = NULL, *scratch = NULL;
double *cpt = NULL;
bool locally_complete = FALSE;
cgdata local_data = { 0 };

  /* allocate a second data table to hold the data that are local to the node's
   * parents. */
  local_data = empty_cgdata(dt.m.nobs, dt.ndcols, dt.ngcols);
  /* allocate the scratch space used to compute the expected value of the
   * distribution of continuous observations. */
  scratch = Calloc1D(dt.m.nobs, sizeof(double));
  /* parent configurations to index the parameters in the probability tables and
   * in the coefficient tables. */
  parcfgs = Calloc1D(dt.m.nobs, sizeof(int));

  /* for each node... */
  for (int i = 0; i < bn.nnodes; i++) {

    /* ... that we want to consider... */
    if (!dt.m.flag[i].fixed)
      continue;

    if (debugging)
      Rprintf("* processing node %s.\n", bn.labels[i]);

    /* ... find out what type of node it is... */
    switch(bn.node_types[i]) {

      /* ... if it's a discrete node... */
      case DNODE:
      case ONODE:

        /* ... get its values from the data and the (conditional) probability
         *  table from the network... */
        dobs = dt.dcol[dt.map[i]];
        cpt = bn.ldists[i].d.cpt;
        locally_complete = dt.m.flag[i].complete;

        if (bn.ldists[i].nparents == 0) {

          /* ... if the node is a root node, the value of the observation is the
           * index in the probability table ... */
          if (locally_complete) {

            for (int j = 0; j < dt.m.nobs; j++)
              loglik[j] += log(cpt[dobs[j] - 1]);

          }/*THEN*/
          else {

            for (int j = 0; j < dt.m.nobs; j++) {

              if (dobs[j] == NA_INTEGER)
                loglik[j] = NA_REAL;
               else
                loglik[j] += log(cpt[dobs[j] - 1]);

            }/*FOR*/

          }/*ELSE*/

        }/*THEN*/
        else {

          /* ... if the node has parents, extract them from the data... */
          cgdata_subset_columns(&dt, &local_data, bn.ldists[i].parents,
             bn.ldists[i].nparents);
          /* ... compute their configurations... */
          c_fast_config(local_data.dcol, local_data.m.nobs, local_data.m.ncols,
              local_data.nlvl, parcfgs, NULL, 0);

          for (int k = 0; k < local_data.m.ncols; k++)
            locally_complete &= local_data.m.flag[k].complete;

          /* ... and the two together index the conditional probability table. */
          if (locally_complete) {

            for (int j = 0; j < dt.m.nobs; j++)
              loglik[j] += log(cpt[CMC(dobs[j] - 1, parcfgs[j],
                               bn.ldists[i].d.dims[0])]);

          }/*THEN*/
          else {

            for (int j = 0; j < dt.m.nobs; j++) {

              if ((dobs[j] == NA_INTEGER) || (parcfgs[j] == NA_INTEGER))
                loglik[j] = NA_REAL;
              else
                loglik[j] += log(cpt[CMC(dobs[j] - 1, parcfgs[j],
                                 bn.ldists[i].d.dims[0])]);

            }/*FOR*/

          }/*ELSE*/

        }/*ELSE*/

        break;

      /* ... if it's a Gaussian node... */
      case GNODE:

        /* ... extract the relevant quantities from the data and the network...*/
        gobs = dt.gcol[dt.map[i]];
        coefs = bn.ldists[i].g.coefs;
        pars = bn.ldists[i].parents;
        sd = bn.ldists[i].g.sd;

        /* ... initialize the log-likelihoods with the intercept... */
        for (int j = 0; j < dt.m.nobs; j++)
          scratch[j] = coefs[0];

        /* ... then add the effects of the parents one by one to compute the
         * expected values... */
        for (int k = 0; k < bn.ldists[i].nparents; k++)
          for (int j = 0; j < dt.m.nobs; j++)
            scratch[j] += dt.gcol[dt.map[pars[k]]][j] * coefs[k + 1];

        /* ... and use them together with the standard error to compute the
         * log-likelihood. */
        for (int j = 0; j < dt.m.nobs; j++)
          loglik[j] += dnorm(gobs[j], scratch[j], sd, TRUE);

        break;

      /* ... if it's a conditional Gaussian node... */
      case CGNODE:

        /* ... extract the relevant quantities from the data and the network...*/
        gobs = dt.gcol[dt.map[i]];
        coefs = bn.ldists[i].cg.coefs;
        ncoefs = bn.ldists[i].cg.ncoefs;
        pars = bn.ldists[i].cg.gparents;
        sds = bn.ldists[i].cg.sd;

        /* ... extract the discrete parents ... */
        cgdata_subset_columns(&dt, &local_data, bn.ldists[i].cg.dparents,
           bn.ldists[i].cg.ndparents);
        /* ... compute their configurations... */
        c_fast_config(local_data.dcol, local_data.m.nobs, local_data.m.ncols,
            local_data.nlvl, parcfgs, NULL, 0);

        /* ... find out if they contain missing values...  */
        locally_complete = TRUE;
        for (int k = 0; k < local_data.m.ncols; k++)
          locally_complete &= local_data.m.flag[k].complete;

        if (locally_complete) {

          /* ... initialize the log-likelihoods with the intercepts... */
          for (int j = 0; j < dt.m.nobs; j++)
            scratch[j] = coefs[CMC(0, parcfgs[j], ncoefs)];

          /* ... then add the effects of the continuous parents one by one to
           * compute the expected values, using the coefficients for the right
           * configuration of the discrete parents ... */
          for (int k = 0; k < bn.ldists[i].cg.ngparents; k++)
            for (int j = 0; j < dt.m.nobs; j++)
              scratch[j] += dt.gcol[dt.map[pars[k]]][j] *
                              coefs[CMC(k + 1, parcfgs[j], ncoefs)];

          /* ... and use them together with the standard error to compute the
           * log-likelihood. */
          for (int j = 0; j < dt.m.nobs; j++)
            loglik[j] += dnorm(gobs[j], scratch[j], sds[parcfgs[j]], TRUE);

        }/*THEN*/
        else {

          /* ... initialize the log-likelihoods to NA if the configuration of
           * the discrete parents is a missing value... */
          for (int j = 0; j < dt.m.nobs; j++) {

            if (parcfgs[j] == NA_INTEGER)
              scratch[j] = NA_REAL;
            else
              scratch[j] = coefs[CMC(0, parcfgs[j], ncoefs)];

          }/*FOR*/

          /* ... but otherwise add the effected of the continuous parents... */
          for (int k = 0; k < bn.ldists[i].cg.ngparents; k++)
            for (int j = 0; j < dt.m.nobs; j++)
              if (!ISNAN(scratch[j]))
                scratch[j] += dt.gcol[dt.map[pars[k]]][j] *
                                coefs[CMC(k + 1, parcfgs[j], ncoefs)];

          /* ... and use them together with the standard error to compute the
           * log-likelihood. */
          for (int j = 0; j < dt.m.nobs; j++) {

            if (ISNAN(scratch[j]))
              loglik[j] = NA_REAL;
            else
              loglik[j] += dnorm(gobs[j], scratch[j], sds[parcfgs[j]], TRUE);

          }/*FOR*/

        }/*ELSE*/

      case ENOFIT:
      default:
        break;

    }/*SWITCH*/

  }/*FOR*/

  Free1D(scratch);
  Free1D(parcfgs);
  FreeCGDT(local_data);

}/*BYSAMPLE_CLGAUSSIAN_LOGLIKELIHOOD*/

/* log-likelihood of a whole sample for a Gaussian network. */
double data_clgaussian_loglikelihood(fitted_bn bn, cgdata dt, double *scratch,
    bool propagate, bool debugging) {

int max_nlvls = 0, cumdim = 0, max_cfgs = 0, ncomplete = 0;
int *pars = NULL, *parcfgs = NULL, ncoefs = 0;
double *gobs = NULL, *coefs = NULL, sd = 0, *sds = NULL;
double loglik = 0, node_loglik = 0;
bool locally_complete = FALSE, early_return = FALSE;
counts1d freq = { 0 };
counts2d freq2 = { 0 };
cgdata local_data = { 0 };

  /* if the data contain missing values for the nodes we are considering, and
   * we propagate them, the log-likelihood is necessarily NA. */
  if (propagate && check_locally_incomplete_data(bn, dt.m, debugging))
    return NA_REAL;

  /* if any of the coefficients of a Gaussian nodes we are considering is NA,
   * or the standard error is NA, then the log-likelihood is NA; this is not
   * necessarily the case for discrete and conditional Gaussian nodes since
   * not all parameters are necessarily involved in the computation. */
  for (int i = 0; i < bn.nnodes; i++) {

    if ((!dt.m.flag[i].fixed) || (bn.node_types[i] != GNODE))
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
          Rprintf("* incomplete data for node %s, the log-likelihood is NA.\n",
            bn.labels[i]);
        return NA_REAL;

      }/*THEN*/

  }/*FOR*/

  /* allocate a second data table to hold the data that are local to the node's
   * parents. */
  local_data = empty_cgdata(dt.m.nobs, dt.ndcols, dt.ngcols);
  /* parent configurations to index the parameters in the probability tables and
   * in the coefficient tables. */
  parcfgs = Calloc1D(dt.m.nobs, sizeof(int));
  /* find out the largest number of levels that a variable can take to bound the
   * size of the contingency tables for the counts. */
  max_nlvls = dt.nlvl[i_which_max(dt.nlvl, dt.ndcols) - 1];
  for (int i = 0; i < bn.nnodes; i++) {

    cumdim = 1;

    switch(bn.node_types[i]) {

      case DNODE:
      case ONODE:

        for (int j = 1; j < bn.ldists[i].d.ndims; j++)
          cumdim *= bn.ldists[i].d.dims[j];
        break;

      case CGNODE:

        cumdim = bn.ldists[i].cg.nconfigs;
        break;

      case ENOFIT:
      default:
        break;

    }/*SWITCH*/

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
    ncomplete = 0;

    /* ... find out what type of node it is... */
    switch(bn.node_types[i]) {

      /* ... if it's a discrete node... */
      case DNODE:
      case ONODE:

        if (bn.ldists[i].nparents == 0) {

          /* ... compute the frequencies... */
          resize_1d_table(dt.nlvl[dt.map[i]], &freq);
          refill_1d_table(dt.dcol[dt.map[i]], &freq, dt.m.nobs);
          /* ... if there are usable locally-complete observations... */
          if (freq.nobs == 0)
            node_loglik = R_NegInf;
          else {

            /* ... combine them with the probabilities to compute the
             * log-likelihood... */
            for (int j = 0; j < freq.llx; j++)
              if (freq.n[j] > 0)
                node_loglik += freq.n[j] * log(bn.ldists[i].d.cpt[j]);
            /* ... and scale the log-likelihood to compensate for any missing
             * values (which will not be propagated as a result). */
            if (freq.nobs < dt.m.nobs)
              node_loglik = node_loglik / freq.nobs * dt.m.nobs;

          }/*ELSE*/

        }/*THEN*/
        else {

          /* ... construct the parent configurations (with a 1-offset to make
           * fill_2d_table() happy)... */
          cgdata_subset_columns(&dt, &local_data, bn.ldists[i].parents,
             bn.ldists[i].nparents);
          c_fast_config(local_data.dcol, local_data.m.nobs, local_data.m.ncols,
              local_data.nlvl, parcfgs, &cumdim, 1);
          /* ... compute the frequencies... */
          resize_2d_table(dt.nlvl[dt.map[i]], cumdim, &freq2);
          refill_2d_table(dt.dcol[dt.map[i]], parcfgs, &freq2, dt.m.nobs);
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
            /* ... and scale the log-likelihood to compensate for any missing
             * values (which will not be propagated as a result). */
            if (freq2.nobs < dt.m.nobs)
              node_loglik = node_loglik / freq2.nobs * dt.m.nobs;

          }/*ELSE*/

        }/*ELSE*/

        break;

      case GNODE:

        /* ... extract the relevant quantities from the data and the network...*/
        gobs = dt.gcol[dt.map[i]];
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
            scratch[j] += dt.gcol[dt.map[pars[k]]][j] * coefs[k + 1];

        /* ... and use them together with the standard error to compute the
         * log-likelihood. */
        for (int j = 0; j < dt.m.nobs; j++) {

          if (ISNAN(gobs[j]) || ISNAN(scratch[j]))
            continue;

          node_loglik += dnorm(gobs[j], scratch[j], sd, TRUE);
          ncomplete++;

        }/*FOR*/

        /* ... and scale it to compensate for any missing values (which will
         * not be propagated as a result), or return -Inf if there are no
         * locally-complete observations. */
        if (ncomplete == 0)
          node_loglik = R_NegInf;
        else if (ncomplete < dt.m.nobs)
          node_loglik = node_loglik / ncomplete * dt.m.nobs;

        break;

      case CGNODE:

        /* ... extract the relevant quantities from the data and the network...*/
        gobs = dt.gcol[dt.map[i]];
        coefs = bn.ldists[i].cg.coefs;
        ncoefs = bn.ldists[i].cg.ncoefs;
        pars = bn.ldists[i].cg.gparents;
        sds = bn.ldists[i].cg.sd;

        /* ... extract the discrete parents ... */
        cgdata_subset_columns(&dt, &local_data, bn.ldists[i].cg.dparents,
           bn.ldists[i].cg.ndparents);
        /* ... compute their configurations... */
        c_fast_config(local_data.dcol, local_data.m.nobs, local_data.m.ncols,
            local_data.nlvl, parcfgs, NULL, 0);

        /* ... find out if they contain missing values...  */
        locally_complete = TRUE;
        for (int k = 0; k < local_data.m.ncols; k++)
          locally_complete &= local_data.m.flag[k].complete;

        if (locally_complete) {

          /* ... initialize the log-likelihoods with the intercepts... */
          for (int j = 0; j < dt.m.nobs; j++)
            scratch[j] = coefs[CMC(0, parcfgs[j], ncoefs)];

          /* ... then add the effects of the continuous parents one by one to
           * compute the expected values, using the coefficients for the right
           * configuration of the discrete parents ... */
          for (int k = 0; k < bn.ldists[i].cg.ngparents; k++)
            for (int j = 0; j < dt.m.nobs; j++)
              scratch[j] += dt.gcol[dt.map[pars[k]]][j] *
                          coefs[CMC(k + 1, parcfgs[j], ncoefs)];

          /* ... and use them together with the standard error to compute the
           * log-likelihood. */
          for (int j = 0; j < dt.m.nobs; j++)
            node_loglik += dnorm(gobs[j], scratch[j], sds[parcfgs[j]], TRUE);

          ncomplete = dt.m.nobs;

        }/*THEN*/
        else {

          /* ... initialize the log-likelihoods to NA if the configuration of
           * the discrete parents is a missing value... */
          for (int j = 0; j < dt.m.nobs; j++) {

            if (parcfgs[j] == NA_INTEGER)
              scratch[j] = NA_REAL;
            else
              scratch[j] = coefs[CMC(0, parcfgs[j], ncoefs)];

          }/*FOR*/

          /* ... but otherwise add the effected of the continuous parents... */
          for (int k = 0; k < bn.ldists[i].cg.ngparents; k++)
            for (int j = 0; j < dt.m.nobs; j++)
              if (!ISNAN(scratch[j]))
                scratch[j] += dt.gcol[dt.map[pars[k]]][j] *
                          coefs[CMC(k + 1, parcfgs[j], ncoefs)];

          /* ... and use them together with the standard error to compute the
           * log-likelihood. */
          for (int j = 0; j < dt.m.nobs; j++) {

            if (ISNAN(gobs[j]) || ISNAN(scratch[j]))
              continue;

            node_loglik += dnorm(gobs[j], scratch[j], sds[parcfgs[j]], TRUE);
            ncomplete++;

          }/*FOR*/

        }/*ELSE*/

        /* ... and scale it to compensate for any missing values (which will
         * not be propagated as a result), or return -Inf if there are no
         * locally-complete observations. */
        if (ncomplete == 0)
          node_loglik = R_NegInf;
        else if (ncomplete < dt.m.nobs)
          node_loglik = node_loglik / ncomplete * dt.m.nobs;

      case ENOFIT:
      default:
        break;

    }/*SWITCH*/

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

  Free1D(parcfgs);
  resize_1d_table(max_nlvls, &freq);
  Free1DTAB(freq);
  resize_2d_table(max_nlvls, max_cfgs, &freq2);
  Free2DTAB(freq2);
  FreeCGDT(local_data);

  return loglik;

}/*DATA_CLGAUSSIAN_LOGLIKELIHOOD*/

