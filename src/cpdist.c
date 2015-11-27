#include "include/rcore.h"
#include "include/sampling.h"
#include "include/dataframe.h"
#include "include/globals.h"

SEXP cpdist_lw(SEXP fitted, SEXP nodes, SEXP n, SEXP fix, SEXP debug) {

int nsims = INT(n);
double *weights = NULL;
SEXP result, simulation, wgt, from;

  /* allocate the scratch space for the simulation. */
  PROTECT(simulation = fit2df(fitted, nsims));
  /* perform the simulation. */
  c_rbn_master(fitted, simulation, n, fix, FALSE);

  if (isTRUE(debug))
    Rprintf("* generated %d samples from the bayesian network.\n", nsims);

  /* compute the weights. */
  PROTECT(wgt = allocVector(REALSXP, nsims));
  weights = REAL(wgt);
  from = getAttrib(fix, R_NamesSymbol);
  c_lw_weights(fitted, simulation, nsims, weights, from, FALSE);

  /* all weights are zero, the event is impossible. */
  if (weights[d_which_max(weights, nsims) - 1] == 0)
    error("all weights are zero, the evidence has probability zero.");

  /* allocate the return value. */
  PROTECT(result = c_dataframe_column(simulation, nodes, FALSE, TRUE));
  minimal_data_frame(result);

  /* set the weights. */
  setAttrib(result, BN_WeightsSymbol, wgt);
  /* allocate and set the method. */
  setAttrib(result, BN_MethodSymbol, mkString("lw"));
  /* allocate and set the class. */
  setAttrib(result, R_ClassSymbol, mkStringVec(2, "bn.cpdist", "data.frame"));

  UNPROTECT(3);

  return result;

}/*CPDIST_LW*/

