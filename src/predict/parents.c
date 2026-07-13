#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/data.table.h"
#include "../core/math.functions.h"
#include "../core/sampling.h"
#include "../fitted/fitted.h"
#include "../math/hyperpoisson.h"
#include "../math/linear.algebra.h"
#include "../scores/scores.h"
#include "predict.h"

/* predict the value of a gaussian node from its (continuous) parents. */
void c_gaussian_predict(ldist ld, double **ccols, int ncont, int nobs,
    double *res, bool debugging) {

double *coefs = ld.g.coefs;

  for (int i = 0; i < nobs; i++) {

    /* compute the mean value for this observation. */
    res[i] = coefs[0];
    for (int j = 0; j < ncont; j++)
      res[i] += ccols[j][i] * coefs[j + 1];

  }/*FOR*/

  if (debugging) {

    if (ncont == 0) {

      Rprintf("  > prediction for all observations is %lf.\n", coefs[0]);

    }/*THEN*/
    else {

      for (int i = 0; i < nobs; i++) {

        Rprintf("  > prediction for observation %d is %lf with predictor:\n",
          i + 1, res[i]);

        Rprintf("    (%lf)", coefs[0]);
        for (int j = 0; j < ncont; j++)
          Rprintf(" + (%lf) * (%lf)", ccols[j][i], coefs[j + 1]);
        Rprintf("\n");

      }/*FOR*/

    }/*ELSE*/

  }/*THEN*/

}/*C_GAUSSIAN_PREDICT*/

/* predict the value of a discrete node from its parent configurations. */
void c_discrete_predict(ldist ld, int *configs, int nconfigs, int nobs,
    int *res, double *pt, bool include_prob, bool debugging) {

int nrow = ld.d.dims[0];
double *cpt = ld.d.cpt;
char **levels = ld.d.levels;
int *iscratch = NULL, *maxima = NULL, *nmax = NULL;
double *dscratch = NULL, *buf = NULL;

  /* allocate the scratch space used by all_max(). */
  iscratch = Calloc1D(nrow, sizeof(int));
  buf = Calloc1D(nrow, sizeof(double));
  dscratch = Calloc1D(nrow * nconfigs, sizeof(double));
  memcpy(dscratch, cpt, nrow * nconfigs * sizeof(double));
  maxima = Calloc1D(nrow * nconfigs, sizeof(int));
  nmax = Calloc1D(nconfigs, sizeof(int));

  /* find out the mode(s) of each configuration. */
  for (int c = 0; c < nconfigs; c++) {

    for (int k = 0; k < nrow; k++)
      iscratch[k] = k + 1;

    nmax[c] = all_max(dscratch + c * nrow, nrow, maxima + c * nrow, iscratch,
                buf);

  }/*FOR*/

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  if (nconfigs == 1) {

    /* node with no parents: a single multinomial distribution shared by all
     * the observations. */
    if (nmax[0] == 1) {

      for (int i = 0; i < nobs; i++)
        res[i] = maxima[0];

      if (debugging) {

        if (res[0] == NA_INTEGER)
          Rprintf("  > prediction for all observations is NA with probabilities:\n");
        else
          Rprintf("  > prediction for all observations is '%s' with probabilities:\n",
            levels[res[0] - 1]);

        Rprintf("  ");
        for (int k = 0; k < nrow; k++)
          Rprintf("  %lf", cpt[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(nobs, nmax[0], res, maxima);

      if (debugging) {

        Rprintf("  > there are %d levels tied for prediction, applying tie breaking.\n",
          nmax[0]);
        Rprintf("  > tied levels are:");
        for (int k = 0; k < nmax[0]; k++)
          Rprintf(" %s", levels[maxima[k] - 1]);
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

    if (include_prob)
      for (int i = 0; i < nobs; i++)
        memcpy(pt + i * nrow, cpt, nrow * sizeof(double));

  }/*THEN*/
  else {

    /* node with parents: one multinomial distribution per configuration. */
    for (int i = 0; i < nobs; i++) {

      if (configs[i] == NA_INTEGER) {

        res[i] = NA_INTEGER;

        if (debugging)
          Rprintf("  > prediction for observation %d is NA because at least one parent is NA.\n", i + 1);

      }/*THEN*/
      else if (nmax[configs[i] - 1] == 0) {

        res[i] = NA_INTEGER;

        if (debugging)
          Rprintf("  > prediction for observation %d is NA because the probabilities are missing.\n", i + 1);

      }/*THEN*/
      else if (nmax[configs[i] - 1] == 1) {

        res[i] = maxima[CMC(0, configs[i] - 1, nrow)];

        if (debugging) {

          if (res[i] == NA_INTEGER)
            Rprintf("  > prediction for observation %d is NA with probabilities:\n", i + 1);
          else
            Rprintf("  > prediction for observation %d is '%s' with probabilities:\n",
              i + 1, levels[res[i] - 1]);

          Rprintf("  ");
          for (int k = 0; k < nrow; k++)
            Rprintf("  %lf", (cpt + nrow * (configs[i] - 1))[k]);
          Rprintf("\n");

        }/*THEN*/

      }/*THEN*/
      else {

        /* break ties: sample with replacement from all the maxima. */
        SampleReplace(1, nmax[configs[i] - 1], res + i,
          maxima + (configs[i] - 1) * nrow);

        if (debugging) {

          Rprintf("  > there are %d levels tied for prediction of observation %d, applying tie breaking.\n",
            nmax[configs[i] - 1], i + 1);
          Rprintf("  > tied levels are:");
          for (int k = 0; k < nmax[configs[i] - 1]; k++)
            Rprintf(" %s", levels[maxima[CMC(k, configs[i] - 1, nrow)] - 1]);
          Rprintf(".\n");

        }/*THEN*/

      }/*ELSE*/

      /* attach the prediction probabilities to the return value. */
      if (include_prob) {

        if (configs[i] == NA_INTEGER)
          for (int k = 0; k < nrow; k++)
            (pt + i * nrow)[k] = NA_REAL;
        else
          memcpy(pt + i * nrow, cpt + nrow * (configs[i] - 1),
            nrow * sizeof(double));

      }/*THEN*/

    }/*FOR*/

  }/*ELSE*/

  /* save the state of the random number generator. */
  PutRNGstate();

  Free1D(iscratch);
  Free1D(buf);
  Free1D(dscratch);
  Free1D(maxima);
  Free1D(nmax);

}/*C_DISCRETE_PREDICT*/

/* predict the value of a conditional gaussian node. */
void c_cgaussian_predict(ldist ld, double **ccols, int ncont, int *configs,
    int nobs, double *res, bool debugging) {

double *beta = ld.cg.coefs, *beta_offset = NULL;

  for (int i = 0; i < nobs; i++) {

    /* is the configuration of the discrete parents defined? */
    if (configs[i] == NA_INTEGER) {

      res[i] = NA_REAL;
      continue;

    }/*THEN*/

    /* find out which conditional regression to use for prediction. */
    beta_offset = beta + (configs[i] - 1) * (ncont + 1);

    /* compute the mean value for this observation. */
    res[i] = beta_offset[0];
    for (int j = 0; j < ncont; j++)
      res[i] += ccols[j][i] * beta_offset[j + 1];

    if (debugging) {

      Rprintf("  > prediction for observation %d is %lf with predictor:\n",
        i + 1, res[i]);

      Rprintf("    (%lf)", beta_offset[0]);
      for (int j = 0; j < ncont; j++)
        Rprintf(" + (%lf) * (%lf)", ccols[j][i], beta_offset[j + 1]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

}/*C_CGAUSSIAN_PREDICT*/

/* predict the (expected) value of a zero-inflated count node. */
void c_zero_inflated_predict(ldist ld, fitted_node_e type, tabular parents, int nobs,
    double *res, bool debugging) {

double *params1 = NULL, *params2 = NULL;

  params1 = Calloc1D(nobs, sizeof(double));
  params2 = Calloc1D(nobs, sizeof(double));

  if (type == ZIHPNODE) {

    double shape = ld.zihp.dispersion;

    /* translate the regression coefficients into canonical parameters. */
    zihp_coefs_to_params(params1, ld.zihp.inflation, params2, ld.zihp.intensity,
      parents);

    for (int i = 0; i < nobs; i++) {

      res[i] = (1 - params1[i]) * hypois_mean(params2[i], shape);

      if (debugging)
        Rprintf("  > prediction for observation %d is %lf with zero-inflation "
          "probability %lf, intensity %lf and dispersion %lf.\n",
          i + 1, res[i], params1[i], params2[i], shape);

    }/*FOR*/

  }/*THEN*/
  else if (type == ZINBNODE) {

    double shape = ld.zinb.failures;

    /* translate the regression coefficients into canonical parameters. */
    zinb_coefs_to_params(params1, ld.zinb.inflation, params2, ld.zinb.prsucc,
      parents);

    for (int i = 0; i < nobs; i++) {

      res[i] = (1 - params1[i]) * shape * params2[i] / (1 - params2[i]);

      if (debugging)
        Rprintf("  > prediction for observation %d is %lf with zero-inflation "
          "probability %lf, success probability %lf, number of failures %lf.\n",
          i + 1, res[i], params1[i], params2[i], shape);

    }/*FOR*/

  }/*THEN*/

  Free1D(params1);
  Free1D(params2);

}/*C_ZIGDAG_PREDICT*/
