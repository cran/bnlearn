#include "include/rcore.h"

SEXP data_frame_finite(SEXP data) {

int i = 0, j = 0, ncols = length(data), nrows = length(VECTOR_ELT(data, 0));
double *xx = 0;
SEXP nodes = getAttrib(data, R_NamesSymbol);

  for (i = 0; i < ncols; i++) {

    xx = REAL(VECTOR_ELT(data, i));

    for (j = 0; j < nrows; j++)
      if (!R_FINITE(xx[j]))
        error("columns %s contains non-finite values.", NODE(i));

  }/*FOR*/

  return R_NilValue;

}/*DATA_FRAME_FINITE*/

SEXP data_type(SEXP data) {

int i = 0, numeric = 0, categorical = 0, ordinal = 0, ncols = length(data);
SEXP column, class, nodes = getAttrib(data, R_NamesSymbol);

  for (i = 0; i  < ncols; i++) {

    column = VECTOR_ELT(data, i);

    switch(TYPEOF(column)) {

      case REALSXP:
        numeric++;
        break;

      case INTSXP:
        class = getAttrib(column, R_ClassSymbol);

        if ((length(class) == 1) &&
            (strcmp(CHAR(STRING_ELT(class, 0)), "factor") == 0))
              categorical++;
        else if ((length(class) == 2) &&
            (strcmp(CHAR(STRING_ELT(class, 0)), "ordered") == 0) &&
            (strcmp(CHAR(STRING_ELT(class, 1)), "factor") == 0))
              ordinal++;
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
const char *cur_node_class = NULL;
SEXP nodes, vars, try_nodes, try_vars, temp, cur_node, cur_node_levels;
SEXP cur_var, cur_var_levels, cur_var_class;

  /* match nodes in the network and variables in the data to the active subset. */
  nodes = getAttrib(fitted, R_NamesSymbol);
  vars = getAttrib(data, R_NamesSymbol);
  PROTECT(try_nodes = match(nodes, subset, 0));
  tn = INTEGER(try_nodes);
  PROTECT(try_vars = match(vars, subset, 0));
  tv = INTEGER(try_vars);

  /* iterate over the variables in the data and cross-check with the nodes. */
  for (i = 0; i < length(subset); i++) {

    cur_var = VECTOR_ELT(data, tv[i] - 1);
    cur_node = VECTOR_ELT(fitted, tn[i] - 1);
    cur_node_class = CHAR(STRING_ELT(getAttrib(cur_node, R_ClassSymbol), 0));

    switch(TYPEOF(cur_var)) {

      case REALSXP:
        if ((strcmp(cur_node_class, "bn.fit.gnode") != 0) &&
            (strcmp(cur_node_class, "bn.fit.cgnode") != 0))
          error("node %s is discrete but variable %s in the data is continuous.",
            NODE(i), NODE(i));
        break;

      case INTSXP:
        if ((strcmp(cur_node_class, "bn.fit.dnode") != 0) &&
            (strcmp(cur_node_class, "bn.fit.onode") != 0))
          error("node %s is continuous but variable %s in the data is discrete.",
            NODE(i), NODE(i));

        /* warn if the variable is categorical and the node is ordinal or vice
         * versa. */
        cur_var_class = getAttrib(cur_var, R_ClassSymbol);

        if ((strcmp(cur_node_class, "bn.fit.dnode") == 0) &&
            (length(cur_var_class) == 2))
          warning("node %s is categorical but variable %s in the data is ordinal.",
            NODE(i), NODE(i));
        else if ((strcmp(cur_node_class, "bn.fit.onode") == 0) &&
                 (length(cur_var_class) == 1))
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
              j + 1, NODE(i), CHAR(STRING_ELT(cur_var_levels, j)),
              CHAR(STRING_ELT(cur_node_levels, j)));

        /* check that levels are the same. */
        break;

      default:
        error("variables must be either numeric, factors or ordered factors.");

    }/*SWITCH*/

  }/*FOR*/

  UNPROTECT(2);

  return R_NilValue;

}/*FITTED_VS_DATA*/

