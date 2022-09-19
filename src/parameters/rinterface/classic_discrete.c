#include "../../include/rcore.h"
#include "../../core/contingency.tables.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/table.h"
#include "../parameters.h"

SEXP classic_discrete_parameters(SEXP data, SEXP node, SEXP parents, SEXP iss,
    SEXP replace_unidentifiable, SEXP missing) {

double *cpt = 0, alpha = 0;
bool replace = isTRUE(replace_unidentifiable);
SEXP nodes_in_order, relevant_data, relevant_df, counts, cptable;

  /* subset the data with the node labels in the right order. */
  PROTECT(nodes_in_order = allocVector(STRSXP, length(parents) + 1));
  SET_STRING_ELT(nodes_in_order, 0, STRING_ELT(node, 0));
  for (int i = 0; i < length(parents); i++)
    SET_STRING_ELT(nodes_in_order, i + 1, STRING_ELT(parents, i));

  PROTECT(relevant_data = c_dataframe_column(data, nodes_in_order, FALSE, TRUE));
  PROTECT(relevant_df = minimal_data_frame(relevant_data));

  /* implement the maximum likelihood estimator as a particular case of the
   * posterior estimator with prior mass equal to zero. */
  if (iss == R_NilValue)
    alpha = 0;
  else
    alpha = NUM(iss);

  /* compute the counts that are the sufficient statistic. */
  PROTECT(counts = minimal_table(relevant_df, missing));

  /* prepare the conditional probability table... */
  PROTECT(cptable = allocVector(REALSXP, length(counts)));
  setAttrib(cptable, R_DimSymbol, getAttrib(counts, R_DimSymbol));
  setAttrib(cptable, R_DimNamesSymbol, getAttrib(counts, R_DimNamesSymbol));
  setAttrib(cptable, R_ClassSymbol, mkString("table"));
  cpt = REAL(cptable);

  /* ... and estimate it. */
  c_classic_discrete_parameters(INTEGER(counts), cpt, nrows(cptable),
      length(cptable) / nrows(cptable), alpha, replace);

  UNPROTECT(5);

  return cptable;

}/*CLASSIC_DISCRETE_PARAMETERS*/
