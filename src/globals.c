
#include "common.h"

SEXP BN_ModelstringSymbol;
SEXP BN_NodesSymbol;
SEXP BN_ProbSymbol;
SEXP BN_ScoreDeltaSymbol;

SEXP c_onLoad() {

  /* initialize symbols in .onLoad(), to do that only once. */
  BN_ModelstringSymbol = install("modelstring");
  BN_NodesSymbol = install("nodes");
  BN_ProbSymbol = install("prob");
  BN_ScoreDeltaSymbol = install("score.delta");

  return R_NilValue;

}/*C_ONLOAD*/
