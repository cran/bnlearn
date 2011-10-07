#include "common.h"

SEXP arcs2elist(SEXP arcs, SEXP nodes, SEXP id, SEXP sublist) {

int i = 0, j = 0, k = 0, nnodes = LENGTH(nodes), narcs = LENGTH(arcs)/2;
int *e = NULL, *coords = NULL, *children = NULL;
int *convert = LOGICAL(id), *sub = LOGICAL(sublist);
SEXP try, elist, edges, temp, temp_name = R_NilValue;

  /* allocate the return value. */
  PROTECT(elist = allocVector(VECSXP, nnodes));
  /* set the node names. */
  setAttrib(elist, R_NamesSymbol, nodes);

  if (*sub > 0) {

    /* allocate and initialize the subset name. */
    PROTECT(temp_name = allocVector(STRSXP, 1));
    SET_STRING_ELT(temp_name, 0, mkChar("edges"));

  }/*THEN*/

  /* allocate the scratch space to keep track of the children of each node. */
  children = alloc1dcont(nnodes);

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  for (i = 0; i < narcs; i++) {

    children[coords[i] - 1]++;

  }/*FOR*/

  for (i = 0; i < nnodes; i++) {

    /* allocate and set up the edge array. */
    if (*convert > 0) {

      PROTECT(edges = allocVector(INTSXP, children[i]));
      e = INTEGER(edges);

    }/*THEN*/
    else {

      PROTECT(edges = allocVector(STRSXP, children[i]));

    }/*ELSE*/

    /* copy the coordinates of the adjacent nodes. */
    for (j = 0, k = 0; j < narcs; j++) {

      if (coords[j] != i + 1)
        continue;

      if (*convert > 0)
        e[k++] = coords[narcs + j];
      else
        SET_STRING_ELT(edges, k++, STRING_ELT(arcs, narcs + j));

      if (k == children[i])
        break;

    }/*FOR*/

    if (*sub > 0) {

      /* allocate and set up the "edge" sublist for graphNEL. */
      PROTECT(temp = allocVector(VECSXP, 1));
      setAttrib(temp, R_NamesSymbol, temp_name);
      SET_VECTOR_ELT(temp, 0, edges);
      SET_VECTOR_ELT(elist, i, temp);
      UNPROTECT(1);

    }/*THEN*/
    else {

      SET_VECTOR_ELT(elist, i, edges);

    }/*ELSE*/

    UNPROTECT(1);

  }/*FOR*/

  if (*sub > 0)
    UNPROTECT(3);
  else
    UNPROTECT(2);

  return elist;

}/*ARCS2ELIST*/

