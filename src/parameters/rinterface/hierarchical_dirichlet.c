#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/table.h"
#include "../../minimal/common.h"
#include "../../math/linear.algebra.h"
#include "../parameters.h"

SEXP hierarchical_dirichlet_parameters(SEXP data, SEXP node, SEXP parents, SEXP
    group, SEXP alpha0, SEXP iss, SEXP missing, SEXP debug) {

int ngroups = 0, cpt_nrows = 0, cpt_ncols = 0;
double *cpt = 0;
long double colsum = 0;
hdstatus err = { 0 };
cmcmap cc = { 0 };
SEXP nodes_in_order, relevant_data, relevant_df, counts, cptable;

  /* subset the data with the node labels in the right order. */
  PROTECT(nodes_in_order = allocVector(STRSXP, length(parents) + 2));
  SET_STRING_ELT(nodes_in_order, 0, STRING_ELT(node, 0));
  for (int i = 0; i < length(parents); i++)
    SET_STRING_ELT(nodes_in_order, i + 1, STRING_ELT(parents, i));
  SET_STRING_ELT(nodes_in_order, length(parents) + 1, STRING_ELT(group, 0));

  PROTECT(relevant_data = c_dataframe_column(data, nodes_in_order, FALSE, TRUE));
  PROTECT(relevant_df = minimal_data_frame(relevant_data));

  /* compute the counts that are the sufficient statistic. */
  PROTECT(counts = minimal_table(relevant_df, missing));

  /* create a contingency table that is convenient to use: one row per
   * configuration of the node and the parents, one column per group. */
  ngroups = NLEVELS(VECTOR_ELT(relevant_df, length(relevant_df) - 1));
  cc.el = INTEGER(counts);
  cc.nrows = length(counts) / ngroups;
  cc.ncols = ngroups;

  /* prepare the conditional probability table... */
  PROTECT(cptable = allocVector(REALSXP, length(counts)));
  setAttrib(cptable, R_DimSymbol, getAttrib(counts, R_DimSymbol));
  setAttrib(cptable, R_DimNamesSymbol, getAttrib(counts, R_DimNamesSymbol));
  setAttrib(cptable, R_ClassSymbol, mkString("table"));
  cpt = REAL(cptable);
  cpt_nrows = nrows(cptable);
  cpt_ncols = length(cptable) / cpt_nrows;

  /* ... estimate it... */
  err = c_hierarchical_dirichlet_parameters(cc, NUM(alpha0), NUM(iss) / cc.ncols,
          isTRUE(debug), cpt);

  /* ... and normalize the columns to sum up to 1. */
  for (int j = 0; j < cpt_ncols; j++) {

    colsum = 0;
    for (int i = 0; i < cpt_nrows; i++)
      colsum += cpt[CMC(i, j, cpt_nrows)];
    for (int i = 0; i < cpt_nrows; i++)
      cpt[CMC(i, j, cpt_nrows)] /= colsum;

  }/*FOR*/

  PrintValue(cptable);

  UNPROTECT(5);

  /* warnings at the end of the function, so that they do not cause leaks even
   * if warnings are transformed into errors via options(). */
  if (err.outer_em_convergence_fail)
    warning("possible convergence failure in the EM outer loop for node %s.",
      CHAR(STRING_ELT(node, 0)));
  if (err.kappa_tau_convergence_fail)
    warning("possible convergence failure in the Newton update for kappa and tau for node %s.",
      CHAR(STRING_ELT(node, 0)));
  if (err.tau_convergence_fail)
    warning("possible convergence failure in the Newton update for tau for node %s.",
      CHAR(STRING_ELT(node, 0)));
  if (err.kappa_convergence_fail)
    warning("possible convergence failure in the Newton update for kappa for node %s.",
      CHAR(STRING_ELT(node, 0)));
  if (err.tau_is_zero)
    warning("tau is zero, restarting the Newton updates for node %s.",
      CHAR(STRING_ELT(node, 0)));

  return cptable;

}/*HIERARCHICAL_DIRICHLET_PARAMETERS*/
