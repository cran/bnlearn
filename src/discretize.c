#include "include/rcore.h"
#include "include/data.table.h"
#include "include/data.frame.h"
#include "include/preprocessing.h"
#include "include/tests.h"

static void interval_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs) {

int i = 0, k = 0;
double min = R_PosInf, max = R_NegInf, delta = 0;

  /* determine the range (maximum and minimum values, assumed finite)...  */
  for (i = 0, min = R_PosInf, max = R_NegInf; i < nobs; i++) {

    if (orig[i] < min)
      min = orig[i];
    if (orig[i] > max)
      max = orig[i];

  }/*FOR*/

  /* ... cut the range in intervals of equal length... */
  delta = (max - min) / nbreaks;

  /* ... derive the factor integer code... */
  for (i = 0; i < nobs; i++) {

    if (orig[i] == min)
      factor[i] = 1;
    else
      factor[i] = ceil((orig[i] - min) / delta);

  }/*FOR*/

  /* ... and save the cutpoints as well. */
  for (k = 0; k < nbreaks; k++)
    cutpoints[k] = min + k * delta;
  cutpoints[nbreaks] = max;

}/*INTERVAL.DISCRETIZATION*/

static void quantile_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs) {

int i = 0, k = 0, lo = 0, hi = 0;
double h = 0, *sorted = NULL;

  /* sort a copy of the data... */
  sorted = Calloc1D(nobs, sizeof(double));
  memcpy(sorted, orig, nobs * sizeof(double));
  R_qsort(sorted, 1, nobs);

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

  /* derive the factor integer code. */
  for (i = 0; i < nobs; i++) {

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

SEXP marginal_discretize(SEXP data, SEXP method, SEXP breaks, SEXP ordered) {

int i = 0, j = 0, max_nbreaks = 0;
int *nbreaks = INTEGER(breaks), *create_ordered = LOGICAL(ordered);
cgdata orig = { 0 };
discretization_e m = discretization_to_enum(CHAR(STRING_ELT(method, 0)));
SEXP discretized, new_factor, new_levels;

  /* store the data in a data table that allows both discrete and continuous
   * variables, with the understanding that continuous variables are not
   * necessarily Gaussian. */
  orig = cgdata_from_SEXP(data, 0, 0);

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > nbreaks[i]) ? max_nbreaks : nbreaks[i];

  if ((m == INTERVAL) || (m == QUANTILE)) {

    /* cutpoints used in the discretization should be stored to produce the
     * labels for the factor levels. */
    double *cutpoints = Calloc1D(max_nbreaks + 1, sizeof(double));

    for (j = 0; j < orig.m.ncols; j++) {

      if (orig.m.flag[j].discrete) {

        /* leave discrete variables alone, just copy them over. */
        SET_VECTOR_ELT(discretized, j, VECTOR_ELT(data, j));
        continue;

      }/*THEN*/

      /* allocate and initialize the factor. */
      PROTECT(new_factor = allocVector(INTSXP, orig.m.nobs));

      if (m == INTERVAL)
        interval_discretization(orig.gcol[orig.map[j]], INTEGER(new_factor),
          nbreaks[j], cutpoints, orig.m.nobs);
      else if (m == QUANTILE)
        quantile_discretization(orig.gcol[orig.map[j]], INTEGER(new_factor),
            nbreaks[j], cutpoints, orig.m.nobs);

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

    Free1D(cutpoints);

  }/*THEN*/

  FreeCGDT(orig, FALSE);

  /* make sure the return value is a data frame. */
  PROTECT(discretized = minimal_data_frame(discretized));

  UNPROTECT(2);

  return discretized;

}/*MARGINAL_DISCRETIZE*/

static void hartemink_discretization(ddata work, int *ibreaks, int *nbreaks,
        int max_nbreaks, double **cutpoints) {

int i = 0;
int ***n = NULL, **ni = NULL, **nj = NULL, *nobs = NULL;

  /* separate marginals for each joint, they will not be identical acrosso
   * contingency tables with missing data. */
  n = Calloc1D(work.m.ncols, sizeof(int **));
  ni = Calloc1D(work.m.ncols, sizeof(int *));
  nj = Calloc1D(work.m.ncols, sizeof(int *));
  nobs = Calloc1D(work.m.ncols, sizeof(int));

  /* for the i-th variable... */
  for (i = 0; i < work.m.ncols; i++) {


  }/*FOR*/

  Free1D(n);
  Free1D(ni);
  Free1D(nj);
  Free1D(nobs);

}/*HARTEMINK_DISCRETIZATION*/

SEXP joint_discretize(SEXP data, SEXP method, SEXP breaks, SEXP ordered,
    SEXP initial_discretization, SEXP initial_breaks) {

int i = 0, j = 0, max_nbreaks = 0;
int *nbreaks = INTEGER(breaks), *ibreaks = INTEGER(initial_breaks);
int *create_ordered = LOGICAL(ordered);
double **all_cutpoints = NULL;
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

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > ibreaks[i]) ? max_nbreaks : ibreaks[i];

  if (m == HARTEMINK) {

    /* store the transformed data after the initial discretization, and pass
     * that to Hartemink's method to keep the whole thing modular. */
    workspace = empty_ddata(orig.m.nobs, orig.m.ncols);
    workspace.m.flag = Calloc1D(orig.m.ncols, sizeof(flags));

    /* cutpoints used in the discretization should be stored to produce the
     * labels for the factor levels. */
    all_cutpoints =
      (double **)Calloc2D(orig.m.ncols, max_nbreaks + 1, sizeof(double));
    all_SEXPs = Calloc1D(orig.m.ncols, sizeof(SEXP));

    for (j = 0; j < orig.m.ncols; j++) {

      /* mark discrete variables as fixed, they should not be modified when
       * levels are collapsed later on. */
      if (orig.m.flag[j].discrete) {

        workspace.col[j] = orig.dcol[orig.map[j]];
        workspace.nlvl[j] = orig.nlvl[orig.map[j]];
        workspace.m.flag[j].fixed = TRUE;
        all_SEXPs[j] = VECTOR_ELT(data, j);

      }/*THEN*/
      else {

        /* perform the initial discretization. */
        PROTECT(new_factor = allocVector(INTSXP, orig.m.nobs));
        workspace.col[j] = INTEGER(new_factor);
        workspace.nlvl[j] = ibreaks[j];

        if (idisc == INTERVAL) {

          interval_discretization(orig.gcol[orig.map[j]], workspace.col[j],
              ibreaks[j], all_cutpoints[j], orig.m.nobs);

        }/*THEN*/
        else if (idisc == QUANTILE) {

          quantile_discretization(orig.gcol[orig.map[j]], workspace.col[j],
              ibreaks[j], all_cutpoints[j], orig.m.nobs);

        }/*THEN*/

        all_SEXPs[j] = new_factor;

      }/*ELSE*/

    }/*FOR*/

    hartemink_discretization(workspace, ibreaks, nbreaks, max_nbreaks,
      all_cutpoints);

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

  FreeCGDT(orig, FALSE);
  Free1D(all_cutpoints);
  Free1D(all_SEXPs);

  /* make sure the return value is a data frame. */
  PROTECT(discretized = minimal_data_frame(discretized));

  UNPROTECT(1 + orig.ngcols);

  return discretized;

}/*JOINT_DISCRETIZE*/
