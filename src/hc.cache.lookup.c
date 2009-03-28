#include "common.h"

#define ARC(i, col) CHAR(STRING_ELT(VECTOR_ELT(cache, col), i))
#define CACHE(i) REAL(VECTOR_ELT(cache, 2))[i]

SEXP cache_lookup(SEXP arc, SEXP cache) {

  int i = 0;
  int nrows = LENGTH(VECTOR_ELT(cache, 0));
  const char *from = CHAR(STRING_ELT(arc, 0));
  const char *to = CHAR(STRING_ELT(arc, 1));
  SEXP res;

  /* allocate and initialize the result. */
  PROTECT(res = allocVector(REALSXP, 1));
  NUM(res) = R_NaN;

  for (i = 0; i < nrows; i++) {

    /* check the first element; if it does not match skip the second one. */
    if (!strcmp(from, ARC(i, 0)) ) {

      /* if the first element matches, check the other one. */
      if (!strcmp(to, ARC(i, 1)) ) {

        /* found cached value, storing it in the result. */
        NUM(res) = CACHE(i);

        /* each arc has at most one entry in the cache, so we can break. */
        break;

      }/*THEN*/

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  return res;

}/*CACHE_LOOKUP*/

