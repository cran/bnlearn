#include "../../include/rcore.h"
#include "../../include/globals.h"
#include "../../core/allocations.h"
#include "../../core/contingency.tables.h"
#include "../../core/sort.h"
#include "../../include/globals.h"
#include "../../core/data.table.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../minimal/common.h"
#include "../../tests/tests.h"
#include "../preprocessing.h"

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

static void copy_cutpoints(double *cutpoints, int nbreaks, SEXP list, int idx) {

SEXP temp;

    PROTECT(temp = allocVector(REALSXP, nbreaks + 1));
    memcpy(REAL(temp), cutpoints, (nbreaks + 1) * sizeof(double));
    SET_VECTOR_ELT(list, idx, temp);
    UNPROTECT(1);

}/*COPY_CUTPOINTS*/

SEXP marginal_discretize(SEXP data, SEXP method, SEXP breaks, SEXP ordered,
    SEXP debug) {

int i = 0, j = 0, max_nbreaks = 0, errcode = 0;
int *nbreaks = INTEGER(breaks), *create_ordered = LOGICAL(ordered);
bool debugging = isTRUE(debug);
const char *cur_node = NULL;
cgdata orig = { 0 };
discretization_e m = discretization_to_enum(CHAR(STRING_ELT(method, 0)));
SEXP discretized, new_factor, new_levels, metadata, complete_nodes;
SEXP attached_cutpoints;

  /* extract the metadata from the data. */
  PROTECT(metadata = getAttrib(data, BN_MetaDataSymbol));
  PROTECT(complete_nodes = getListElement(metadata, "complete.nodes"));

  /* store the data in a data table that allows both discrete and continuous
   * variables, with the understanding that continuous variables are not
   * necessarily Gaussian. */
  orig = cgdata_from_SEXP(data, 0, 0);
  meta_copy_names(&(orig.m), 0, data);
  meta_init_flags(&(orig.m), 0, complete_nodes, R_NilValue);

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > nbreaks[i]) ? max_nbreaks : nbreaks[i];

  /* set up the cutpoints, which will be attached to the return value. */
  PROTECT(attached_cutpoints = allocVector(VECSXP, orig.m.ncols));
  setAttrib(attached_cutpoints, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

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
            nbreaks[j], cutpoints, orig.m.nobs, orig.m.flag[j].complete,
            debugging);

      }/*THEN*/

      /* check whether the discretization was successful, and produce an error
       * instead of returning a data frame if not. */
      if (errcode) {

        Free1D(cutpoints);
        FreeCGDT(orig);
        UNPROTECT(4);

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

      /* save the cutpoints. */
      copy_cutpoints(cutpoints, nbreaks[j], attached_cutpoints, j);

      /* save the factor into the return value. */
      SET_VECTOR_ELT(discretized, j, new_factor);
      UNPROTECT(2);

    }/*FOR*/

  }/*THEN*/

  Free1D(cutpoints);
  FreeCGDT(orig);

  /* make sure the return value is a data frame. */
  PROTECT(discretized = minimal_data_frame(discretized));
  /* attach the cutpoints to make the discretization reproducible. */
  setAttrib(discretized, BN_CutpointsSymbol, attached_cutpoints);

  UNPROTECT(5);

  return discretized;

}/*MARGINAL_DISCRETIZE*/

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
SEXP discretized, new_factor, new_levels, metadata, complete_nodes;
SEXP *all_SEXPs = NULL;
SEXP attached_cutpoints;

  /* extract the metadata from the data. */
  PROTECT(metadata = getAttrib(data, BN_MetaDataSymbol));
  PROTECT(complete_nodes = getListElement(metadata, "complete.nodes"));

  /* store the data in a data table that allows both discrete and continuous
   * variables, with the understanding that continuous variables are not
   * necessarily Gaussian. */
  orig = cgdata_from_SEXP(data, 0, 0);
  meta_copy_names(&(orig.m), 0, data);
  meta_init_flags(&(orig.m), 0, complete_nodes, R_NilValue);

  /* set up the return value. */
  PROTECT(discretized = allocVector(VECSXP, orig.m.ncols));
  setAttrib(discretized, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

  for (i = 0; i < orig.m.ncols; i++)
    max_nbreaks = (max_nbreaks > ibreaks[i]) ? max_nbreaks : ibreaks[i];

  /* set up the cutpoints, which will be attached to the return value. */
  PROTECT(attached_cutpoints = allocVector(VECSXP, orig.m.ncols));
  setAttrib(attached_cutpoints, R_NamesSymbol, getAttrib(data, R_NamesSymbol));

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
            ibreaks[j], all_cutpoints[j], orig.m.nobs, orig.m.flag[j].complete,
            debugging);

      }/*THEN*/

      /* check whether the discretization was successful, try again with fewer
       * cutpoints if possible. */
      if (errcode) {

        if (ibreaks[j] == nbreaks[j]) {

          Free2D(all_cutpoints, orig.m.ncols);
          Free1D(all_SEXPs);
          FreeDDT(workspace);
          FreeCGDT(orig);
          UNPROTECT(3);

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

        /* save the cutpoints. */
        copy_cutpoints(all_cutpoints[j], nbreaks[j], attached_cutpoints, j);

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
  /* attach the cutpoints to make the discretization reproducible. */
  setAttrib(discretized, BN_CutpointsSymbol, attached_cutpoints);

  UNPROTECT(5 + orig.ngcols);

  return discretized;

}/*JOINT_DISCRETIZE*/
