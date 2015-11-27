#include "include/rcore.h"

SEXP BN_ModelstringSymbol;
SEXP BN_NodesSymbol;
SEXP BN_ProbSymbol;
SEXP BN_MethodSymbol;
SEXP BN_WeightsSymbol;
SEXP TRUESEXP, FALSESEXP;

SEXP c_onLoad() {

  /* initialize symbols in .onLoad(), to do that only once. */
  BN_ModelstringSymbol = install("modelstring");
  BN_NodesSymbol = install("nodes");
  BN_ProbSymbol = install("prob");
  BN_MethodSymbol = install("method");
  BN_WeightsSymbol = install("weights");
  TRUESEXP = ScalarLogical(TRUE);
  R_PreserveObject(TRUESEXP);
  FALSESEXP = ScalarLogical(FALSE);
  R_PreserveObject(FALSESEXP);

  return R_NilValue;

}/*C_ONLOAD*/

SEXP c_onUnload() {

  R_ReleaseObject(TRUESEXP);
  R_ReleaseObject(FALSESEXP);

  return R_NilValue;

}/*C_ONUNLOAD*/

