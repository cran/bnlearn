#include <R.h>
#include <Rinternals.h>

#define ARC(i,col) CHAR(STRING_ELT(arcs, i + col * nrows))
#define NODE(i) CHAR(STRING_ELT(nodes, i))
#define AMAT(i,j) INTEGER(amat)[i + j * nrows]

SEXP arcs2amat(SEXP arcs, SEXP nodes) {

  int i = 0, j = 0, k = 0;
  int nrows = LENGTH(arcs) / 2;
  int dims = LENGTH(nodes);
  SEXP res, dimnames;

  /* allocate the adjacency matrix. */
  PROTECT(res = allocMatrix(INTSXP, dims, dims));

  /* allocate rownames and colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, nodes);
  SET_VECTOR_ELT(dimnames, 1, nodes);
  setAttrib(res, R_DimNamesSymbol, dimnames);

  /* initialize the adjacency matrix. */
  memset(INTEGER(res), '\0', sizeof(int) * dims * dims);

  /* nothing to do if there are no arcs. */
  if (nrows == 0) {

    UNPROTECT(2);
    return res;

  }/*THEN*/

  /* iterate over the arcs. */
  for (k = 0; k < nrows; k++) {

    /* iterate over labels to match the parent's. */
    for (i = 0; i < dims; i++) {

      if (!strcmp(NODE(i), ARC(k, 0)) ) {

        /* iterate over labels to match the parent's. */
        for (j = 0; j < dims; j++) {

          if (!strcmp(NODE(j), ARC(k, 1)) ) {

            INTEGER(res)[i + j * dims] = 1;

          }/*THEN*/

        }/*FOR*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

  return res;

}/*ARCS2AMAT*/

SEXP amat2arcs(SEXP amat, SEXP nodes) {

  int i = 0, j = 0, k = 0;
  int nrows = LENGTH(nodes);
  int narcs = 0;
  SEXP arcs, dimnames, colnames;

  /* count the number of arcs in the adjacency matrix. */
  for (i = 0; i < nrows; i++) {

    for (j = 0; j < nrows; j++) {

      if (AMAT(i, j) == 1) narcs++;

    }/*FOR*/

  }/*FOR*/

  /* allocate colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* if there are no arcs, return an empty arc list. */
  if (narcs == 0) {

    /* allocate an empty arc list. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc list. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

  }/*ELSE*/

  /* fill the arc list from the adjacency matrix. */
  for (i = 0; i < nrows; i++) {

    for (j = 0; j < nrows; j++) {

      /* colnames and rownames are completely ignored. This kills some corner
           cases present in the old R code.  */
      if (AMAT(i, j) == 1) {

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

