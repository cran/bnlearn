#include "common.h"

#define ARC(i,col) CHAR(STRING_ELT(arcs, i + col * nrows))
#define NODE(i) CHAR(STRING_ELT(nodes, i))

SEXP arcs2amat(SEXP arcs, SEXP nodes) {

  int i = 0, j = 0, k = 0;
  int nrows = LENGTH(arcs) / 2;
  int dims = LENGTH(nodes);
  int *res;
  SEXP result, dimnames;

  /* allocate the adjacency matrix. */
  PROTECT(result = allocMatrix(INTSXP, dims, dims));
  res = INTEGER(result);

  /* allocate rownames and colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, nodes);
  SET_VECTOR_ELT(dimnames, 1, nodes);
  setAttrib(result, R_DimNamesSymbol, dimnames);

  /* initialize the adjacency matrix. */
  memset(res, '\0', sizeof(int) * dims * dims);

  /* nothing to do if there are no arcs. */
  if (nrows == 0) {

    UNPROTECT(2);
    return result;

  }/*THEN*/

  /* iterate over the arcs. */
  for (k = 0; k < nrows; k++) {

    /* iterate over labels to match the parent's. */
    for (i = 0; i < dims; i++) {

      if (!strcmp(NODE(i), ARC(k, 0)) ) {

        /* iterate over labels to match the parent's. */
        for (j = 0; j < dims; j++) {

          if (!strcmp(NODE(j), ARC(k, 1)) ) {

            res[i + j * dims] = 1;

          }/*THEN*/

        }/*FOR*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

  return result;

}/*ARCS2AMAT*/

SEXP amat2arcs(SEXP amat, SEXP nodes) {

  int i = 0, j = 0, k = 0;
  int nrows = LENGTH(nodes), narcs = 0;
  int *a = INTEGER(amat);
  SEXP arcs, dimnames, colnames;

  /* count the number of arcs in the adjacency matrix. */
  for (i = 0; i < nrows; i++) {

    for (j = 0; j < nrows; j++) {

      if (a[CMC(i, j, nrows)] == 1) narcs++;

    }/*FOR*/

  }/*FOR*/

  /* allocate colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    /* allocate an empty arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

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

  UNPROTECT(3);

  return arcs;

}/*AMAT2ARCS*/

