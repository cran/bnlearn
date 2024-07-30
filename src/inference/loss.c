#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/data.frame.h"
#include "../minimal/common.h"
#include "../include/globals.h"
#include "../core/sets.h"
#include "../fitted/fitted.h"
#include "../math/linear.algebra.h"

/* classification error of a single node as a loss function. */
SEXP class_err(SEXP reference, SEXP predicted) {

int i = 0, dropped = 0, ndata = length(reference);
int *r = INTEGER(reference), *p = INTEGER(predicted);
double err = 0;

  /* count how many elements differ (this assumes the levels of the two factors
   * are the same and in the same order); NAs are dropped. */
  for (i = 0; i < ndata; i++) {

    if ((r[i] == NA_INTEGER) || (p[i] == NA_INTEGER))
      dropped++;
    else if (r[i] != p[i])
      err++;

  }/*FOR*/

  /* rescale into a probability. */
  if (ndata > dropped)
    err /= (ndata - dropped);
  else
    err = NA_REAL;

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because of missing values.", dropped);

  return ScalarReal(err);

}/*CLASS_ERR*/
