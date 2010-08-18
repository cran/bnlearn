#include "common.h"

/* convert an arc set to an adjacency matrix. */
SEXP arcs2amat(SEXP arcs, SEXP nodes) {

int k = 0, nrows = LENGTH(arcs) / 2, dims = LENGTH(nodes);
int *res = NULL, *coords = NULL;
SEXP result, dimnames, try;

  /* allocate and initialize the adjacency matrix. */
  PROTECT(result = allocMatrix(INTSXP, dims, dims));
  res = INTEGER(result);
  memset(res, '\0', sizeof(int) * dims * dims);

  /* allocate rownames and colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, nodes);
  SET_VECTOR_ELT(dimnames, 1, nodes);
  setAttrib(result, R_DimNamesSymbol, dimnames);

  /* nothing to do if there are no arcs. */
  if (nrows == 0) {

    UNPROTECT(2);
    return result;

  }/*THEN*/

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  /* iterate over the arcs. */
  for (k = 0; k < nrows; k++)
    res[CMC(coords[k] - 1, coords[k + nrows] - 1, dims)] = 1;

  UNPROTECT(3);

  return result;

}/*ARCS2AMAT*/

/* convert an adjacency matrix to an arc set. */
SEXP amat2arcs(SEXP amat, SEXP nodes) {

int i = 0, j = 0, k = 0, nrows = LENGTH(nodes), narcs = 0;
int *a = INTEGER(amat);
SEXP arcs;

  /* count the number of arcs in the adjacency matrix. */
  for (i = 0; i < nrows; i++) {

    for (j = 0; j < nrows; j++) {

      if (a[CMC(i, j, nrows)] == 1) narcs++;

    }/*FOR*/

  }/*FOR*/

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    /* allocate an empty arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    finalize_arcs(arcs);

    UNPROTECT(1);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    finalize_arcs(arcs);

  }/*ELSE*/

  /* fill the arc set from the adjacency matrix. */
  for (i = 0; i < nrows; i++) {

    for (j = 0; j < nrows; j++) {

      /* colnames and rownames are completely ignored. This kills some corner
           cases present in the old R code.  */
      if (a[CMC(i, j, nrows)] == 1) {

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

