#include "../../include/rcore.h"
#include "../../core/data.table.h"
#include "../../fitted/fitted.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/strings.h"
#include "../predict.h"

/* predict new observations from a naive or tree-augmented naive bayes model. */
SEXP naivepred(SEXP fitted, SEXP data, SEXP training, SEXP prior, SEXP prob,
    SEXP debug) {

int tr_nlvl = 0, tr_id = INT(training) - 1;
int *res = NULL, *map = NULL;
double *pr = NULL, *pt = NULL;
bool debugging = isTRUE(debug), include_prob = isTRUE(prob);
SEXP result, data_nodes, fitted_nodes, probtab = R_NilValue, try, tr_levels;
fitted_bn bn = fitted_network_from_SEXP(fitted);
tabular dt = tabular_from_SEXP(data, 0, 0);

  /* cache the node labels. */
  PROTECT(data_nodes = getAttrib(data, R_NamesSymbol));
  PROTECT(fitted_nodes = getAttrib(fitted, R_NamesSymbol));
  PROTECT(try = match(data_nodes, fitted_nodes, 0));
  map = INTEGER(try);
  for (int i = 0; i < length(try); i++)
    map[i]--;

  /* get the training variable and its levels. */
  tr_nlvl = bn.ldists[tr_id].d.dims[0];
  /* get the prior distribution. */
  pr = REAL(prior);

  if (debugging) {

    Rprintf("* the prior distribution for the target variable is:\n");
    PrintValue(prior);

  }/*THEN*/

  /* allocate the return value. */
  PROTECT(result = allocVector(INTSXP, dt.m.nobs));
  res = INTEGER(result);

  /* allocate and initialize the table of the posterior probabilities. */
  if (include_prob) {

    PROTECT(probtab = allocMatrix(REALSXP, tr_nlvl, dt.m.nobs));
    pt = REAL(probtab);
    memset(pt, '\0', dt.m.nobs * tr_nlvl * sizeof(double));

  }/*THEN*/

  /* run the prediction on the bnlearn data structures. */
  c_naivepred(dt, bn, tr_id, tr_nlvl, pr, map, res, pt, include_prob,
    debugging);

  /* add back the attributes and the class to the return value. */
  PROTECT(tr_levels = allocVector(STRSXP, tr_nlvl));
  for (int k = 0; k < tr_nlvl; k++)
    SET_STRING_ELT(tr_levels, k, mkChar(bn.ldists[tr_id].d.levels[k]));
  setAttrib(result, R_LevelsSymbol, tr_levels);

  switch(bn.node_types[tr_id]) {

    case DNODE:
      setAttrib(result, R_ClassSymbol, mkString("factor"));
      break;

    case ONODE:
      setAttrib(result, R_ClassSymbol, mkStringVec(2, "ordered", "factor"));
      break;

    default:
      error("unknown node type.");

  }/*SWITCH*/

  if (include_prob > 0) {

    /* set the levels of the target variable as rownames. */
    setDimNames(probtab, tr_levels, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(6);

  }/*THEN*/
  else {

    UNPROTECT(5);

  }/*ELSE*/

  FreeTAB(dt);
  FreeFittedBN(bn);

  return result;

}/*NAIVEPRED*/

