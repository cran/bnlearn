#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/contingency.tables.h"
#include "../core/data.table.h"
#include "../core/sort.h"
#include "../include/globals.h"
#include "../tests/tests.h"
#include "preprocessing.h"

int interval_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs, bool debugging) {

int i = 0, k = 0;
double min = R_PosInf, max = R_NegInf, delta = 0;

  if (debugging)
    Rprintf("  > discretizing in %d levels.\n", nbreaks);

  /* determine the range (maximum and minimum values, assumed finite)...  */
  for (i = 0, min = R_PosInf, max = R_NegInf; i < nobs; i++) {

    if (orig[i] < min)
      min = orig[i];
    if (orig[i] > max)
      max = orig[i];

  }/*FOR*/

  /* ... cut the range in intervals of equal length... */
  delta = (max - min) / nbreaks;

  if (debugging)
    Rprintf("  > the range is [%lf, %lf], the interval length is %lf.\n",
        min, max, delta);

  /* ... derive the factor integer code (missing values are handled
   * appropriately because ceil(NA) is still NA)... */
  for (i = 0; i < nobs; i++) {

    if (orig[i] == min)
      factor[i] = 1;
    else
      factor[i] = (int)ceil((orig[i] - min) / delta);

  }/*FOR*/

  /* ... and save the cutpoints as well. */
  for (k = 0; k < nbreaks; k++)
    cutpoints[k] = min + k * delta;
  cutpoints[nbreaks] = max;

  /* set the return status to an error if cutpoints are not distinct, that is,
   * if they produce zero-length intervals. */
  for (k = 1; k < nbreaks; k++)
    if (fabs(cutpoints[k] - cutpoints[k - 1]) < MACHINE_TOL)
      return -1;

  return 0;

}/*INTERVAL.DISCRETIZATION*/

int quantile_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs, bool complete, bool debugging) {

int i = 0, k = 0, nonNA = 0, lo = 0, hi = 0;
double h = 0, *sorted = NULL;

  if (debugging)
    Rprintf("  > discretizing in %d levels.\n", nbreaks);

  /* sort a copy of the data, disregarding missing values... */
  sorted = Calloc1D(nobs, sizeof(double));

  if (complete) {

    memcpy(sorted, orig, nobs * sizeof(double));
    nonNA = nobs;

  }/*THEN*/
  else {

    for (i = 0; i < nobs; i++)
      if (!ISNAN(orig[i]))
        sorted[nonNA++] = orig[i];

  }/*ELSE*/

  d_sort(sorted, NULL, nonNA);

  if (debugging)
    Rprintf("  > the range is [%lf, %lf].\n", sorted[0], sorted[nonNA - 1]);

  /* ... compute the indexes of the cutpoints... */
  for (k = 0; k < nbreaks; k++)
    cutpoints[k] = (nonNA - 1) * ((double)k / nbreaks);
  cutpoints[k] = nonNA - 1;

  /* ... compute the cutpoints using quantile() type 7. */
  for (k = 1; k < nbreaks; k++) {

    lo = (int)floor(cutpoints[k]);
    hi = (int)ceil(cutpoints[k]);

    if ((cutpoints[k]) > lo && (sorted[hi] != sorted[lo])) {

      h = (cutpoints[k] - lo);
      cutpoints[k] = (1 - h) * sorted[lo] + h * sorted[hi];

    }/*THEN*/
    else {

      cutpoints[k] = sorted[lo];

    }/*ELSE*/

  }/*FOR*/
  cutpoints[0] = sorted[0];
  cutpoints[nbreaks] = sorted[nobs - 1];

  /* set the return status to an error if cutpoints are not distinct, that is,
   * if they produce zero-length intervals. */
  for (k = 1; k < nbreaks; k++) {

    if (fabs(cutpoints[k] - cutpoints[k - 1]) < MACHINE_TOL) {

      Free1D(sorted);
      return -1;

    }/*THEN*/

  }/*FOR*/

  /* derive the factor integer code. */
  for (i = 0; i < nobs; i++) {

    /* preserve missing values. */
    if (ISNAN(orig[i])) {

      factor[i] = NA_INTEGER;
      continue;

    }/*THEN*/

    for (k = nbreaks - 1; k >= 0; k--) {

      if (orig[i] > cutpoints[k]) {

        factor[i] = k + 1;
        break;

      }/*THEN*/

      /* this is to catch the minimum. */
      factor[i] = 1;

    }/*FOR*/

  }/*FOR*/

  Free1D(sorted);

  return 0;

}/*QUANTILE_DISCRETIZATION*/

void hartemink_discretization(tabular work, int *nbreaks, double **cutpoints,
    bool debugging) {

int i = 0, j = 0, k = 0, max_nlvl = 0, n = work.m.nobs, index_best = 0;
long double cumulated = 0, candidate = 0, current_best = 0;
bool all_done = FALSE;
counts2d *counts = NULL;

  /* contingecy tables become smaller in each iteration: allocate them once with
   * their initial (maximum) dimensions and then gradually make them smaller,
   * zeroing them before each use. */
  counts = Calloc1D(work.m.ncols, sizeof(counts2d));
  for (i = 0; i < work.m.ncols; i++)
    max_nlvl = (max_nlvl < work.nlvl[i]) ? work.nlvl[i] : max_nlvl;
  for (i = 0; i < work.m.ncols; i++)
    counts[i] = new_2d_table(max_nlvl, max_nlvl, TRUE);

  do {

    /* for the i-th variable... */
    for (i = 0, all_done = TRUE; i < work.m.ncols; i++) {

      /* ... if the variable was not already discrete in the first place... */
      if (work.m.flag[i].fixed) {

        if (debugging)
          Rprintf("* skipping variable %s.\n", work.m.names[i]);

        continue;

      }/*THEN*/

      if (work.nlvl[i] > nbreaks[i])
        all_done = FALSE;
      else
        continue;

      if (debugging)
        Rprintf("* Hartemink discretization of variable %s.\n", work.m.names[i]);

      /* ... and the j-th variable... */
      for (j = 0, cumulated = 0; j < work.m.ncols; j++) {

        if (i == j)
          continue;

        /* ... size and fill the contingency tables... */
        resize_2d_table(work.nlvl[i], work.nlvl[j], &(counts[j]));
        refill_2d_table(work.dcol[i], work.dcol[j], &(counts[j]), n);

        /* ... tally up the mutual informations. */
        cumulated += mi_kernel(counts[j]) / counts[j].nobs;

      }/*FOR*/

      if (debugging)
        Rprintf("  > mutual information is %lf.\n", (double)cumulated);

      /* for each level of the i-th variable... */
      for (k = 0, current_best = 0, index_best = 0; k < work.nlvl[i] - 1; k++) {

        if (debugging)
          Rprintf("  > collapsing [%g, %g] and [%g, %g], ",
            cutpoints[i][k], cutpoints[i][k + 1],
            cutpoints[i][k + 1], cutpoints[i][k + 2]);

        /* ... tally up the mutual informations after collapsing the level and the
         * next one... */
        for (j = 0, candidate = 0; j < work.m.ncols; j++)
          if (i != j)
            candidate += mi_kernel_collapsed(counts[j], k) / counts[j].nobs;

        if (debugging)
          Rprintf("mutual information is now %lf.\n", (double)candidate);

        /* ... and pick the level which increases the mutual information the
         * least. */
        if (candidate > current_best) {

          current_best = candidate;
          index_best = k;

        }/*THEN*/

      }/*FOR*/

      if (debugging)
        Rprintf("  @ best collapse is [%g, %g] and [%g, %g] with mutual information %lf.\n",
          cutpoints[i][index_best], cutpoints[i][index_best + 1],
          cutpoints[i][index_best + 1], cutpoints[i][index_best + 2],
          (double)current_best);

      /* remove the cutpoint in between the now-merged levels (remember that
       * there is one more cutpoint than the number of breaks. */
      for (k = index_best + 1; k < work.nlvl[i]; k++)
        cutpoints[i][k] = cutpoints[i][k + 1];
      cutpoints[i][work.nlvl[i]] = NA_REAL;

      /* adjust the integer coding of the factors in the data frame. */
      for (int l = 0; l < work.m.nobs; l++)
        if (work.dcol[i][l] > index_best + 1)
          work.dcol[i][l]--;
      work.nlvl[i]--;

    }/*FOR*/

  } while (!all_done);

  for (i = 0; i < work.m.ncols; i++) {

    resize_2d_table(max_nlvl, max_nlvl, &(counts[i]));
    Free2DTAB(counts[i]);

  }/*FOR*/
  Free1D(counts);

}/*HARTEMINK_DISCRETIZATION*/

