#include <R.h>
#include <Rinternals.h>

/* utility macros */

#define isTRUE(logical) LOGICAL(logical)[0] == TRUE

/* utility functions  */

SEXP getListElement(SEXP list, char *str);

/* from arcs2amat.c  */

SEXP arcs2amat(SEXP arcs, SEXP nodes);

