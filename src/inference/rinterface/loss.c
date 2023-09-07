#include "../../include/rcore.h"
#include "../../minimal/strings.h"
#include "../loss.h"

/* compute the entropy loss, for use in cross-validation. */
SEXP entropy_loss(SEXP fitted, SEXP data, SEXP debug) {

int ndata = length(VECTOR_ELT(data, 0));
double loss = 0, effective = 0;
SEXP loss_info = R_NilValue, keep = R_NilValue;

  /* always use all the variables in computing the entropy. */
  PROTECT(keep = getAttrib(fitted, R_NamesSymbol));

  loss = c_entropy_loss(fitted, data, ndata, FALSE, NULL, &effective, keep,
                  TRUE, TRUE, isTRUE(debug));

  PROTECT(loss_info = allocVector(VECSXP, 2));
  setAttrib(loss_info, R_NamesSymbol, mkStringVec(2, "loss",  "effective.size"));
  SET_VECTOR_ELT(loss_info, 0, ScalarReal(loss));
  SET_VECTOR_ELT(loss_info, 1, ScalarReal(effective));

  UNPROTECT(2);

  return loss_info;

}/*ENTROPY_LOSS*/

