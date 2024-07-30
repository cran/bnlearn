#include "../../include/rcore.h"
#include "../../include/sampling.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../include/globals.h"
#include "../../core/math.functions.h"

/* sample observations with the associated likelihood weights. */
SEXP cpdist_lw(SEXP fitted, SEXP nodes, SEXP n, SEXP fix, SEXP debug) {

int nsims = INT(n), max_id = 0;
double *weights = NULL;
SEXP result, simulation, wgt, from;

  /* allocate the scratch space for the simulation. */
  PROTECT(simulation = fit2df(fitted, nsims));
  /* generate the random observations. */
  c_rbn_master(fitted, simulation, n, fix, TRUE, FALSE);

  if (isTRUE(debug))
    Rprintf("* generated %d samples from the bayesian network.\n", nsims);

  /* compute the weights. */
  PROTECT(wgt = allocVector(REALSXP, nsims));
  weights = REAL(wgt);
  PROTECT(from = getAttrib(fix, R_NamesSymbol));
  c_lw_weights(fitted, simulation, nsims, weights, from, FALSE);

  /* if all weights are zero or NA, the evidence is making it impossible to
   * generate a set of valid random observations. */
  max_id = d_which_max(weights, nsims);

  if (max_id == NA_INTEGER)
    error("all weights are NA, the probability of the evidence is impossible to compute.");
  if (weights[d_which_max(weights, nsims) - 1] == 0)
    error("all weights are zero, the evidence has probability zero.");

  /* prepare the return value with all the attributes. */
  PROTECT(result = c_dataframe_column(simulation, nodes, FALSE, TRUE));
  minimal_data_frame(result);
  setAttrib(result, BN_WeightsSymbol, wgt);
  setAttrib(result, BN_MethodSymbol, mkString("lw"));
  setAttrib(result, R_ClassSymbol, mkStringVec(2, "bn.cpdist", "data.frame"));

  UNPROTECT(4);

  return result;

}/*CPDIST_LW*/

