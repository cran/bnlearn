#include "common.h"

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
SEXP column, class;

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
          error("variables must be either numeric, factors or ordered factors.");

        break;

      default:
        error("variables must be either numeric, factors or ordered factors.");

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
