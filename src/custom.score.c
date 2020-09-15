#include "include/rcore.h"
#include "include/data.frame.h"

/* evaluate a user-provided, custom decomposable score function. */
double custom_score_function(SEXP target, SEXP x, SEXP data, SEXP custom_fn,
    SEXP custom_args, bool debugging) {

SEXP nodes_info, target_info, parents, call, args_iterator, result;

  /* allocate and populate the pairlist to be valuated. */
  PROTECT(args_iterator = call = allocList(5));
  SET_TYPEOF(call, LANGSXP);
  /* first slot, the function name. */
  SETCAR(args_iterator, custom_fn);
  args_iterator = CDR(args_iterator);
  /* second slot, the label of the target node. */
  SETCAR(args_iterator, target);
  args_iterator = CDR(args_iterator);
  /* third slot, the labels of the parents. */
  nodes_info = getListElement(x, "nodes");
  target_info = getListElement(nodes_info, (char *)CHAR(STRING_ELT(target, 0)));
  parents = getListElement(target_info, "parents");
  SETCAR(args_iterator, parents);
  args_iterator = CDR(args_iterator);
  /* fourth slot, the data. */
  SETCAR(args_iterator, data);
  args_iterator = CDR(args_iterator);
  /* fifth slot, the optional arguments passed as a list. */
  SETCAR(args_iterator, custom_args);
  /* evaluate the custom score function. */
  PROTECT(result = eval(call, R_GlobalEnv));

  /* the return value must be a scalar, real number. */
  if ((TYPEOF(result) != REALSXP) || (length(result) != 1))
    error("the score for node %s must be a scalar, real value.",
      CHAR(STRING_ELT(target, 0)));

  UNPROTECT(2);

  return NUM(result);

}/*CUSTOM_SCORE_FUNCTION*/
