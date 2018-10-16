#include "include/rcore.h"
#include "include/fitted.h"

SEXP data_frame_finite(SEXP data) {

int i = 0, j = 0, ncol = length(data), nrow = length(VECTOR_ELT(data, 0));
double *xx = 0;
SEXP nodes = getAttrib(data, R_NamesSymbol);

  for (i = 0; i < ncol; i++) {

    xx = REAL(VECTOR_ELT(data, i));

    for (j = 0; j < nrow; j++)
      if (!R_FINITE(xx[j]) && !ISNAN(xx[j]))
        error("columns %s contains non-finite values.", NODE(i));

  }/*FOR*/

  return R_NilValue;

}/*DATA_FRAME_FINITE*/

SEXP count_observed_values(SEXP data) {

int i = 0, j = 0, ncol = length(data), nrow = length(VECTOR_ELT(data, 0));
int *rr = NULL, *cc = NULL, *temp_integer = NULL;
double *temp_real = NULL;
SEXP counts, rows, cols, temp;

  PROTECT(counts = allocVector(VECSXP, 2));
  setAttrib(counts, R_NamesSymbol, mkStringVec(2, "rows", "columns"));
  PROTECT(rows = allocVector(INTSXP, nrow));
  PROTECT(cols = allocVector(INTSXP, ncol));
  setAttrib(cols, R_NamesSymbol, getAttrib(data, R_NamesSymbol));
  SET_VECTOR_ELT(counts, 0, rows);
  SET_VECTOR_ELT(counts, 1, cols);
  rr = INTEGER(rows);
  cc = INTEGER(cols);
  memset(rr, '\0', nrow * sizeof(int));
  memset(cc, '\0', ncol * sizeof(int));

  for (j = 0; j < ncol; j++) {

    temp = VECTOR_ELT(data, j);

    switch(TYPEOF(temp)) {

       case REALSXP:

         temp_real = REAL(temp);
         for (i = 0; i < nrow; i++) {

           rr[i] += !ISNAN(temp_real[i]);
           cc[j] += !ISNAN(temp_real[i]);

         }/*FOR*/

         break;

       case INTSXP:

         temp_integer = INTEGER(temp);
         for (i = 0; i < nrow; i++) {

           rr[i] += (temp_integer[i] != NA_INTEGER);
           cc[j] += (temp_integer[i] != NA_INTEGER);

         }/*FOR*/

         break;

    }/*SWITCH*/

  }/*FOR*/

  UNPROTECT(3);

  return counts;

}/*COUNT_OBSERVED_VALUES*/

SEXP data_type(SEXP data) {

int i = 0, numeric = 0, categorical = 0, ordinal = 0, ncol = length(data);
SEXP column, nodes = getAttrib(data, R_NamesSymbol);

  for (i = 0; i  < ncol; i++) {

    column = VECTOR_ELT(data, i);

    switch(TYPEOF(column)) {

      case REALSXP:
        numeric++;
        break;

      case INTSXP:
        if (c_is(column, "ordered"))
          ordinal++;
        else if (c_is(column, "factor"))
          categorical++;
        else
          error("variable %s is not supported in bnlearn (type: %s).",
            NODE(i), type2char(TYPEOF(column)));

        break;

      default:
          error("variable %s is not supported in bnlearn (type: %s).",
            NODE(i), type2char(TYPEOF(column)));

    }/*SWITCH*/

  }/*FOR*/

  if (numeric > 0) {

    if ((categorical == 0) && (ordinal == 0))
      return mkString("continuous");
    else
      return mkString("mixed-cg");

  }/*THEN*/
  else {

    if ((categorical == 0) && (ordinal > 0))
      return mkString("ordered");
    else if ((categorical > 0) && (ordinal == 0))
      return mkString("factor");
    else
      return mkString("mixed-do");

  }/*ELSE*/

}/*DATA_TYPE*/

SEXP fitted_vs_data (SEXP fitted, SEXP data, SEXP subset) {

int i = 0, j = 0, *tn = NULL, *tv = NULL;
fitted_node_e cur_node_type = ENOFIT;
SEXP nodes, vars, try_nodes, try_vars, temp, cur_node, cur_node_levels;
SEXP cur_var, cur_var_levels, cur_var_class;

  /* match nodes in the network and variables in the data to the active subset. */
  PROTECT(nodes = getAttrib(fitted, R_NamesSymbol));
  PROTECT(vars = getAttrib(data, R_NamesSymbol));
  PROTECT(try_nodes = match(nodes, subset, 0));
  tn = INTEGER(try_nodes);
  PROTECT(try_vars = match(vars, subset, 0));
  tv = INTEGER(try_vars);

  /* iterate over the variables in the data and cross-check with the nodes. */
  for (i = 0; i < length(subset); i++) {

    cur_var = VECTOR_ELT(data, tv[i] - 1);
    cur_node = VECTOR_ELT(fitted, tn[i] - 1);
    cur_node_type = r_fitted_node_label(cur_node);

    switch(TYPEOF(cur_var)) {

      case REALSXP:
        if ((cur_node_type != GNODE) && (cur_node_type != CGNODE))
          error("node %s is discrete but variable %s in the data is continuous.",
            NODE(i), NODE(i));
        break;

      case INTSXP:
        if ((cur_node_type != DNODE) && (cur_node_type != ONODE))
          error("node %s is continuous but variable %s in the data is discrete.",
            NODE(i), NODE(i));

        /* warn if the variable is categorical and the node is ordinal or vice
         * versa. */
        cur_var_class = getAttrib(cur_var, R_ClassSymbol);

        if ((cur_node_type == DNODE) && (length(cur_var_class) == 2))
          warning("node %s is categorical but variable %s in the data is ordinal.",
            NODE(i), NODE(i));
        else if ((cur_node_type == ONODE) && (length(cur_var_class) == 1))
          warning("node %s is ordinal but variable %s in the data is categorical.",
            NODE(i), NODE(i));

        /* check that levels are the same. */
        cur_var_levels = getAttrib(cur_var, R_LevelsSymbol);
        temp = getAttrib(getListElement(cur_node, "prob"), R_DimNamesSymbol);
        cur_node_levels = VECTOR_ELT(temp, 0);

        if (length(cur_node_levels) != length(cur_var_levels))
          error("'%s' has different number of levels in the node and in the data.",
            NODE(i));

        /* check that levels are the same. */
        for (j = 0; j < length(cur_node_levels); j++)
          if (strcmp(CHAR(STRING_ELT(cur_var_levels, j)),
                     CHAR(STRING_ELT(cur_node_levels, j))) != 0)
            error("level %d of %s is '%s' in the node and '%s' in the data.",
              j + 1, NODE(i), CHAR(STRING_ELT(cur_node_levels, j)),
              CHAR(STRING_ELT(cur_var_levels, j)));

        break;

      default:
        error("variables must be either numeric, factors or ordered factors.");

    }/*SWITCH*/

  }/*FOR*/

  UNPROTECT(4);

  return R_NilValue;

}/*FITTED_VS_DATA*/

