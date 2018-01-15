#include "include/rcore.h"
#include "include/matrix.h"
#include "include/bn.h"

/* convert an arc set to an adjacency matrix. */
SEXP arcs2amat(SEXP arcs, SEXP nodes) {

int k = 0, nrow = length(arcs) / 2, dims = length(nodes);
int *res = NULL, *coords = NULL;
SEXP result, try;

  /* allocate and initialize the adjacency matrix. */
  PROTECT(result = allocMatrix(INTSXP, dims, dims));
  res = INTEGER(result);
  memset(res, '\0', sizeof(int) * dims * dims);

  /* set rownames and colnames to the node labels. */
  setDimNames(result, nodes, nodes);

  /* nothing to do if there are no arcs. */
  if (nrow == 0) {

    UNPROTECT(1);
    return result;

  }/*THEN*/

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  /* iterate over the arcs. */
  for (k = 0; k < nrow; k++)
    res[CMC(coords[k] - 1, coords[k + nrow] - 1, dims)] = 1;

  UNPROTECT(2);

  return result;

}/*ARCS2AMAT*/

/* convert an adjacency matrix to an arc set. */
SEXP amat2arcs(SEXP amat, SEXP nodes) {

int i = 0, j = 0, k = 0, nrow = length(nodes), narcs = 0;
int *a = INTEGER(amat);
SEXP arcs;

  /* count the number of arcs in the adjacency matrix. */
  for (i = 0; i < nrow; i++) {

    for (j = 0; j < nrow; j++) {

      if (a[CMC(i, j, nrow)] == 1) narcs++;

    }/*FOR*/

  }/*FOR*/

  /* allocate the arc set and set the column names. */
  PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
  setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    UNPROTECT(1);
    return arcs;

  }/*THEN*/

  /* fill the arc set from the adjacency matrix. */
  for (i = 0; i < nrow; i++) {

    for (j = 0; j < nrow; j++) {

      /* colnames and rownames are completely ignored to avoid hitting some
       * corner cases present in the old R code.  */
      if (a[CMC(i, j, nrow)] == 1) {

         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, i));
         SET_STRING_ELT(arcs, k + 1 * narcs, STRING_ELT(nodes, j));
         k++;

      }/*THEN*/

      /* no more arcs, get out of both loops. */
      if (k == narcs) goto end;

    }/*FOR*/

  }/*FOR*/

end:

  UNPROTECT(1);

  return arcs;

}/*AMAT2ARCS*/

