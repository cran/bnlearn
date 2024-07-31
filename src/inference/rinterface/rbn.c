#include "../../include/rcore.h"
#include "../../include/sampling.h"
#include "../../minimal/data.frame.h"

/* generate random observations from a bayesian network. */
SEXP rbn_master(SEXP fitted, SEXP n, SEXP fix, SEXP add_metadata, SEXP debug) {

SEXP result;

  /* allocate the return value. */
  PROTECT(result = fit2df(fitted, INT(n)));
  /* generate the random observations. */
  c_rbn_master(fitted, result, n, fix, isTRUE(add_metadata), isTRUE(debug));

  UNPROTECT(1);

  return result;

}/*RBN_MASTER*/

