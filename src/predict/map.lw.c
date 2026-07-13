#include "../include/rcore.h"
#include "../core/data.table.h"
#include "../core/math.functions.h"
#include "../fitted/fitted.h"
#include "../include/sampling.h"
#include "../math/linear.algebra.h"
#include "../sanitization/data.h"
#include "predict.h"

static double posterior_mean(double *x, double *wgt, int n, int *drop,
    bool debugging) {

int k = 0;
long double wsum = 0, wtot = 0;

  for (k = 0; k < n; k++) {

    /* c_rbn_master() may generate NAs, c_lw_weights() can generate NA and NaNs
       as well, disregard and print a warning. */
    if (ISNAN(x[k]) || ISNAN(wgt[k])) {

      (*drop)++;

    }/*THEN*/
    else {

      wsum += x[k] * wgt[k];
      wtot += wgt[k];

    }/*ELSE*/

  }/*FOR*/

  /* if all weights are zero, the predicted value is NA.*/
  if (wtot == 0)
    wsum = NA_REAL;
  else
    wsum /= wtot;

  if (debugging) {

    if (ISNAN(wsum))
      Rprintf("  > prediction is NA.\n");
    else
      Rprintf("  > prediction is %lf.\n", (double)wsum);

  }/*THEN*/

  return (double)wsum;

}/*POSTERIOR_MEAN*/

static int posterior_mode(int *x, double *wgt, int n, long double *counts,
    char **levels, int nlvls, int *drop, bool debugging) {

int k = 0, max_prob = 0;

  memset(counts, '\0', nlvls * sizeof(long double));

  for (k = 0; k < n; k++) {

    /* c_rbn_master() may generate NAs, and c_lw_weights() can generate NaNs as
     * well, disregard and print a warning. */
    if ((x[k] == NA_INTEGER) || (ISNAN(wgt[k])))
      (*drop)++;
    else
      counts[x[k] - 1] += wgt[k];

  }/*FOR*/

  max_prob = ld_which_max(counts, nlvls);

  /* if all weights are zero, the predicted value is NA.*/
  if (counts[max_prob - 1] == 0)
    max_prob = NA_INTEGER;

  if (debugging) {

    Rprintf("  > prediction is '%s' with weight sums:\n",
      max_prob == NA_INTEGER ? "NA" : levels[max_prob - 1]);
    for (k = 0; k < nlvls; k++)
      Rprintf("%lf ", (double)(counts[k]));
    Rprintf("\n");

  }/*THEN*/

  return max_prob;

}/*POSTERIOR_MODE*/

/* predict the values of one or more variables given one or more variables by
 * likelihood weighting. */
void c_mappred(tabular dt, fitted_bn bn, int node_id, bool target_discrete,
    int nlvls, int nev, bool *ev_discrete, void **varptrs, void **evptrs,
    double *wgt, int nsims, long double *lvls_counts, void *pred, void *res,
    double *pt, bool include_prob, int *drop, bool debugging, SEXP fitted,
    SEXP cpdist, SEXP from, fixed_node *fixed) {

int i = 0, j = 0;
long double lvls_tot = 0;
tabular cpdist_tab = tabular_from_SEXP(cpdist, 0, 0);

  /* iterate over the observations. */
  for (i = 0; i < dt.m.nobs; i++) {

    /* copy the values of the current observation into the evidence: evptrs[]
     * point into the fixed_node structure that conditions the generation. */
    for (j = 0; j < nev; j++) {

      if (ev_discrete[j])
        *((int *)evptrs[j]) = ((int *)varptrs[j])[i];
      else
        *((double *)evptrs[j]) = ((double *)varptrs[j])[i];

    }/*FOR*/

    if (debugging) {

      Rprintf("* predicting observation %d conditional on:\n", i + 1);
      for (j = 0; j < nev; j++) {

        if (ev_discrete[j])
          Rprintf("  %s = %d", CHAR(STRING_ELT(from, j)), *((int *)evptrs[j]));
        else
          Rprintf("  %s = %lf", CHAR(STRING_ELT(from, j)),
            *((double *)evptrs[j]));

      }/*FOR*/
      Rprintf("\n");

    }/*THEN*/

    /* generate samples from the conditional posterior distributions. */
    c_rbn_master(bn, cpdist_tab, fixed, FALSE);
    /* attach the metadata that c_lw_weights() relies on. */
    add_simulation_metadata(cpdist, nsims);
    /* compute the weights. */
    tabular lw_dt = lw_data_from_SEXP(cpdist, fitted, from);
    c_lw_weights(bn, lw_dt, wgt, FALSE);
    FreeTAB(lw_dt);

    /* compute the posterior estimate. */
    if (target_discrete) {

      /* pick the most frequent value. */
      ((int *)res)[i] = posterior_mode((int *)pred, wgt, nsims, lvls_counts,
                          bn.ldists[node_id].d.levels, nlvls, drop, debugging);

      /* compute the posterior probabilities on the right scale, to attach
       * them to the return value. */
      if (include_prob) {

        for (j = 0, lvls_tot = 0; j < nlvls; j++) {

          pt[CMC(j, i, nlvls)] = lvls_counts[j];
          lvls_tot += lvls_counts[j];

        }/*FOR*/

        for (j = 0; j < nlvls; j++)
          pt[CMC(j, i, nlvls)] /= lvls_tot;

      }/*THEN*/

    }/*THEN*/
    else {

      /* average the predicted values. */
      ((double *)res)[i] = posterior_mean((double *)pred, wgt, nsims,
                             drop, debugging);

    }/*ELSE*/

  }/*FOR*/

  FreeTAB(cpdist_tab);

}/*C_MAPPRED*/
