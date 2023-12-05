#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/contingency.tables.h"
#include "../include/globals.h"
#include "../core/data.table.h"
#include "../minimal/data.frame.h"
#include "../minimal/strings.h"
#include "../tests/tests.h"
#include "preprocessing.h"

static int interval_discretization(double *orig, int *factor, int nbreaks,
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

static int quantile_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs, bool debugging) {

int i = 0, k = 0, lo = 0, hi = 0;
double h = 0, *sorted = NULL;

  if (debugging)
    Rprintf("  > discretizing in %d levels.\n", nbreaks);

  /* sort a copy of the data... */
  sorted = Calloc1D(nobs, sizeof(double));
  memcpy(sorted, orig, nobs * sizeof(double));
  R_qsort(sorted, 1, nobs);

  if (debugging)
    Rprintf("  > the range is [%lf, %lf].\n", sorted[0], sorted[nobs - 1]);

  /* ... compute the indexes of the cutpoints... */
  for (k = 0; k < nbreaks; k++)
    cutpoints[k] = (nobs - 1) * ((double)k / nbreaks);
  cutpoints[k] = nobs - 1;

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

/* transform an array of cutpoints into label strings that resemble those
 * produced by cut() in R. */
static SEXP cutpoints_to_levels(double *cutpoints, int nbreaks) {

int i = 0;
char buf[100];
SEXP levels;

  PROTECT(levels = allocVector(STRSXP, nbreaks));

  for (i = 0; i < nbreaks; i++) {

    snprintf(buf, 100, "%s%g,%g]", (i == 0) ? "[" : "(",
      cutpoints[i], cutpoints[i + 1]);

    SET_STRING_ELT(levels, i, mkChar(buf));

  }/*FOR*/

  UNPROTECT(1);

  return levels;

}/*CUTPOINTS_TO_LEVELS*/

SEXP marginal_discretize(SEXP data, SEXP method, SEXP breaks, SEXP ordered,
    SEXP debug) {

int i = 0, j = 0, max_nbreaks = 0, errcode = 0;
int *nbreaks = INTEGER(breaks), *create_ordered = LOGICAL(ordered);
bool debugging = isTRUE(debug);
const char *cur_node = NULL;
cgdata orig = { 0 };
discretization_e m = discretization_to_enum(CHAR(STRING_ELT(method, 0)));
SEXP discretized, new_factor, new_levels;

  /* store the data in a data table that allows both discrete and continuous
   * variables, with the understanding that continuous variables are not
   * necessarily Gaussian. */
  orig = cgdata_from_SEXP(data, 0, 0);
  meta_copy_names(&(orig.m), 0, data);

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > nbreaks[i]) ? max_nbreaks : nbreaks[i];

  /* cutpoints used in the discretization should be stored to produce the
   * labels for the factor levels. */
  double *cutpoints = Calloc1D(max_nbreaks + 1, sizeof(double));

  if ((m == INTERVAL) || (m == QUANTILE)) {

    for (j = 0; j < orig.m.ncols; j++) {

      cur_node = orig.m.names[j];

      if (debugging)
        Rprintf("* %s discretization of variable %s.\n",
          m == INTERVAL ? "interval" : "quantile", cur_node);

      if (orig.m.flag[j].discrete) {

        /* leave discrete variables alone, just copy them over. */
        SET_VECTOR_ELT(discretized, j, VECTOR_ELT(data, j));

        if (debugging)
          Rprintf("  > skipping variable %s, already discrete.\n", cur_node);

        continue;

      }/*THEN*/

      /* allocate and initialize the factor. */
      PROTECT(new_factor = allocVector(INTSXP, orig.m.nobs));

      if (m == INTERVAL) {

        errcode =
          interval_discretization(orig.gcol[orig.map[j]], INTEGER(new_factor),
            nbreaks[j], cutpoints, orig.m.nobs, debugging);

      }/*THEN*/
      else if (m == QUANTILE) {

        errcode =
          quantile_discretization(orig.gcol[orig.map[j]], INTEGER(new_factor),
            nbreaks[j], cutpoints, orig.m.nobs, debugging);

      }/*THEN*/

      /* check whether the discretization was successful, and produce an error
       * instead of returning a data frame if not. */
      if (errcode) {

        Free1D(cutpoints);
        FreeCGDT(orig);
        UNPROTECT(3);

        error("discretizing variable %s into %d levels produced zero-length intervals.",
          cur_node, nbreaks[j]);

      }/*THEN*/

      /* set the levels and the class label. */
      PROTECT(new_levels = cutpoints_to_levels(cutpoints, nbreaks[j]));
      setAttrib(new_factor, R_LevelsSymbol, new_levels);
      if (create_ordered[j])
        setAttrib(new_factor, R_ClassSymbol, mkStringVec(2, "ordered", "factor"));
      else
        setAttrib(new_factor, R_ClassSymbol, mkString("factor"));

      /* save the factor into the return value. */
      SET_VECTOR_ELT(discretized, j, new_factor);
      UNPROTECT(2);

    }/*FOR*/


  }/*THEN*/

  Free1D(cutpoints);
  FreeCGDT(orig);

  /* make sure the return value is a data frame. */
  PROTECT(discretized = minimal_data_frame(discretized));

  UNPROTECT(2);

  return discretized;

}/*MARGINAL_DISCRETIZE*/

static void hartemink_discretization(ddata work, int *nbreaks,
    double **cutpoints, bool debugging) {

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
        refill_2d_table(work.col[i], work.col[j], &(counts[j]), n);

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
        if (work.col[i][l] > index_best + 1)
          work.col[i][l]--;
      work.nlvl[i]--;

    }/*FOR*/

  } while (!all_done);

  for (i = 0; i < work.m.ncols; i++) {

    resize_2d_table(max_nlvl, max_nlvl, &(counts[i]));
    Free2DTAB(counts[i]);

  }/*FOR*/
  Free1D(counts);

}/*HARTEMINK_DISCRETIZATION*/

SEXP joint_discretize(SEXP data, SEXP method, SEXP breaks, SEXP ordered,
    SEXP initial_discretization, SEXP initial_breaks, SEXP debug) {

int i = 0, j = 0, max_nbreaks = 0, errcode = 0;
int *nbreaks = INTEGER(breaks), *ibreaks = INTEGER(initial_breaks);
int *create_ordered = LOGICAL(ordered);
double **all_cutpoints = NULL;
const char *cur_node = NULL;
bool debugging = isTRUE(debug);
ddata workspace = { 0 };
cgdata orig = { 0 };
discretization_e m = discretization_to_enum(CHAR(STRING_ELT(method, 0)));
discretization_e idisc =
  discretization_to_enum(CHAR(STRING_ELT(initial_discretization, 0)));
SEXP discretized, new_factor, new_levels;
SEXP *all_SEXPs = NULL;

  /* store the data in a data table that allows both discrete and continuous
   * variables, with the understanding that continuous variables are not
   * necessarily Gaussian. */
  orig = cgdata_from_SEXP(data, 0, 0);
  meta_copy_names(&(orig.m), 0, data);

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > ibreaks[i]) ? max_nbreaks : ibreaks[i];

  if (m == HARTEMINK) {

    /* store the transformed data after the initial discretization, and pass
     * that to Hartemink's method to keep the whole thing modular. */
    workspace = empty_ddata(orig.m.nobs, orig.m.ncols);
    meta_copy_names(&(workspace.m), 0, data);

    /* cutpoints used in the discretization should be stored to produce the
     * labels for the factor levels. */
    all_cutpoints =
      (double **)Calloc2D(orig.m.ncols, max_nbreaks + 1, sizeof(double));
    all_SEXPs = Calloc1D(orig.m.ncols, sizeof(SEXP));

    for (j = 0; j < orig.m.ncols; j++) {

      cur_node = orig.m.names[orig.map[j]];

      /* mark discrete variables as fixed, they should not be modified when
       * levels are collapsed later on. */
      if (orig.m.flag[j].discrete) {

        workspace.col[j] = orig.dcol[orig.map[j]];
        workspace.nlvl[j] = orig.nlvl[orig.map[j]];
        workspace.m.flag[j].fixed = TRUE;
        all_SEXPs[j] = VECTOR_ELT(data, j);

        continue;

      }/*THEN*/

      if (debugging)
        Rprintf("* %s discretization of variable %s.\n",
          idisc == INTERVAL ? "interval" : "quantile", cur_node);

      /* perform the initial discretization using intervals or quantiles. */
      PROTECT(new_factor = allocVector(INTSXP, orig.m.nobs));
      workspace.col[j] = INTEGER(new_factor);

try_again:

      workspace.nlvl[j] = ibreaks[j];

      if (idisc == INTERVAL) {

        errcode =
          interval_discretization(orig.gcol[orig.map[j]], workspace.col[j],
            ibreaks[j], all_cutpoints[j], orig.m.nobs, debugging);

      }/*THEN*/
      else if (idisc == QUANTILE) {

        errcode =
          quantile_discretization(orig.gcol[orig.map[j]], workspace.col[j],
            ibreaks[j], all_cutpoints[j], orig.m.nobs, debugging);

      }/*THEN*/

      /* check whether the discretization was successful, try again with fewer
       * cutpoints if possible. */
      if (errcode) {

        if (ibreaks[j] == nbreaks[j]) {

          Free2D(all_cutpoints, orig.m.ncols);
          Free1D(all_SEXPs);
          FreeDDT(workspace);
          FreeCGDT(orig);
          UNPROTECT(2);

          error("discretizing node %s produces zero-length intervals even with %d cutpoints.",
            cur_node, nbreaks[j]);

        }/*THEN*/
        else {

          if (debugging)
            Rprintf("  > reducing cutpoints from %d to %d.\n",
              ibreaks[j], ibreaks[j] - 1);

          ibreaks[j]--;
          goto try_again;

        }/*ELSE*/

      }/*THEN*/

      all_SEXPs[j] = new_factor;

    }/*FOR*/

    /* run Hartemink's algorithm to collapse the intervals produced by interval
     * or quantile discretization above. */
    hartemink_discretization(workspace, nbreaks, all_cutpoints, debugging);

    /* construct the labels of the levels and set them. */
    for (j = 0; j < workspace.m.ncols; j++) {

      if (orig.m.flag[j].gaussian) {

        PROTECT(new_levels = cutpoints_to_levels(all_cutpoints[j], nbreaks[j]));
        setAttrib(all_SEXPs[j], R_LevelsSymbol, new_levels);
        if (create_ordered[j])
          setAttrib(all_SEXPs[j], R_ClassSymbol, mkStringVec(2, "ordered", "factor"));
        else
          setAttrib(all_SEXPs[j], R_ClassSymbol, mkString("factor"));

        UNPROTECT(1);

      }/*THEN*/

      SET_VECTOR_ELT(discretized, j, all_SEXPs[j]);

    }/*FOR*/

  }/*THEN*/
  else {

    error("unknown joint discretization method.");

  }/*ELSE*/

  FreeCGDT(orig);
  FreeDDT(workspace);
  Free2D(all_cutpoints, orig.m.ncols);
  Free1D(all_SEXPs);

  /* make sure the return value is a data frame. */
  PROTECT(discretized = minimal_data_frame(discretized));

  UNPROTECT(2 + orig.ngcols);

  return discretized;

}/*JOINT_DISCRETIZE*/
