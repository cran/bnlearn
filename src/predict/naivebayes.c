#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/data.table.h"
#include "../core/math.functions.h"
#include "../core/sampling.h"
#include "../fitted/fitted.h"
#include "../math/linear.algebra.h"

/* predict new observations from a naive or tree-augmented naive bayes model. */
void c_naivepred(tabular dt, fitted_bn bn, int tr_id, int tr_nlvl, double *pr,
    int *map, int *res, double *pt, bool include_prob, bool debugging) {

int nmax = 0, cur_nlvl = 0, idx = 0;
int *iscratch = NULL, *maxima = NULL;
double *cpt = NULL, *scratch = NULL, *buf = NULL;
bool incomplete_observation = FALSE;

  /* allocate the scratch space used to compute posterior probabilities. */
  scratch = Calloc1D(tr_nlvl, sizeof(double));
  buf = Calloc1D(tr_nlvl, sizeof(double));

  /* create the vector of indexes. */
  iscratch = Calloc1D(tr_nlvl, sizeof(int));

  /* allocate the array for the indexes of the maxima. */
  maxima = Calloc1D(tr_nlvl, sizeof(int));

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  /* for each observation... */
  for (int i = 0; i < dt.m.nobs; i++) {

    if (debugging)
      Rprintf("* predicting the value of observation %d.\n", i + 1);

    /* ... check whether it is a complete observation... */
    incomplete_observation = FALSE;

    for (int j = 0; j < dt.m.ncols; j++) {

      if (dt.dcol[j][i] == NA_INTEGER) {

        /* (... mark the observation...) */
        incomplete_observation = TRUE;
        /* (... nuke the class probabilities...) */
        for (int k = 0; k < tr_nlvl; k++)
          scratch[k] = NA_REAL;
        /* (... and skip all processing.)*/
        goto wrap_up;

      }/*THEN*/

    }/*FOR*/

    /* ... reset the scratch space and the indexes array... */
    for (int k = 0; k < tr_nlvl; k++) {

      scratch[k] = log(pr[k]);
      iscratch[k] = k + 1;

    }/*FOR*/

    /* ... and for each conditional probability table... */
    for (int j = 0; j < bn.nnodes; j++) {

      /* ... skip the training variable... */
      if (j == tr_id)
        continue;

      cpt = bn.ldists[j].d.cpt;
      cur_nlvl = bn.ldists[j].d.dims[0];

      /* ... (this is the root node of the Chow-Liu tree, or any node in a naive
       * Bayes model) ... */
      if (bn.ldists[j].nparents == 1) {

        /* ... and for each row of the conditional probability table... */
        for (int k = 0; k < tr_nlvl; k++) {

          if (debugging) {

            Rprintf("  > node %s: picking cell %d (%d, %d) from the CPT (p = %lf).\n",
              bn.labels[j], CMC(dt.dcol[map[j]][i] - 1, k, cur_nlvl),
              dt.dcol[map[j]][i], k + 1,
              cpt[CMC(dt.dcol[map[j]][i] - 1, k, cur_nlvl)]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] +=
            log(cpt[CMC(dt.dcol[map[j]][i] - 1, k, cur_nlvl)]);

        }/*FOR*/

      }/*THEN*/
      else {

        int pr_id = (bn.ldists[j].parents[0] == tr_id) ?
                     bn.ldists[j].parents[1] : bn.ldists[j].parents[0];

        /* ... and for each row of the conditional probability table... */
        for (int k = 0; k < tr_nlvl; k++) {

          /* (the first dimension corresponds to the current node [X], the second
           * to the training node [Y], the third to the only parent of the current
           * node [Z]; CMC coordinates are computed as X + Y * NX + Z * NX * NY. */
          idx = (dt.dcol[map[j]][i] - 1) + k * cur_nlvl +
                  (dt.dcol[map[pr_id]][i] - 1) * cur_nlvl * tr_nlvl;

          if (debugging) {

            Rprintf("  > node %s: picking cell %d (%d, %d, %d) from the CPT (p = %lf).\n",
              bn.labels[j], idx, dt.dcol[map[j]][i], k + 1,
              dt.dcol[map[pr_id]][i], cpt[idx]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] += log(cpt[idx]);

        }/*FOR*/

      }/*ELSE*/

    }/*FOR*/

    /* find out the mode(s). */
    nmax = all_max(scratch, tr_nlvl, maxima, iscratch, buf);

wrap_up:

    /* compute the posterior probabilities on the right scale, to attach them
     * to the return value. */
    if (include_prob) {

      /* copy the log-probabilities from scratch. */
      memcpy(pt + i * tr_nlvl, scratch, tr_nlvl * sizeof(double));

      if (!incomplete_observation) {

        double sum = 0;

        /* transform log-probabilities into plain probabilities. */
        for (int k = 0; k < tr_nlvl; k++)
          pt[i * tr_nlvl + k] =
            exp(pt[i * tr_nlvl + k] - scratch[maxima[0] - 1]);
        for (int k = 0; k < tr_nlvl; k++)
          sum += pt[i * tr_nlvl + k];

        /* rescale them to sum up to 1. */
        for (int k = 0; k < tr_nlvl; k++)
          pt[i * tr_nlvl + k] /= sum;

      }/*THEN*/

    }/*THEN*/

    if (incomplete_observation) {

      res[i] = NA_INTEGER;

      if (debugging)
        Rprintf(" > prediction for observation %d is NA because at least one predictor is NA.\n", i + 1);

    }/*THEN*/
    else if (nmax == 0) {

      res[i] = NA_INTEGER;

      if (debugging)
        Rprintf("  > prediction for observation %d is NA because the probabilities are missing.\n", i + 1);

    }/*THEN*/
    else if (nmax == 1) {

      res[i] = maxima[0];

      if (debugging) {

        Rprintf("  @ prediction for observation %d is '%s' with (log-)posterior:\n",
          i + 1, bn.ldists[tr_id].d.levels[res[i] - 1]);

        Rprintf("  ");
        for (int k = 0; k < tr_nlvl; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(1, nmax, res + i, maxima);

      if (debugging) {

        Rprintf("  @ there are %d levels tied for prediction of observation %d, applying tie breaking.\n", nmax, i + 1);

        Rprintf("  ");
        for (int k = 0; k < tr_nlvl; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

        Rprintf("  @ tied levels are:");
        for (int k = 0; k < nmax; k++)
          Rprintf(" %s", bn.ldists[tr_id].d.levels[maxima[k] - 1]);
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  /* save the state of the random number generator. */
  PutRNGstate();

  Free1D(scratch);
  Free1D(buf);
  Free1D(iscratch);
  Free1D(maxima);

}/*C_NAIVEPRED*/

