#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/common.h"
#include "../minimal/strings.h"

SEXP arcs2ualist(SEXP arcs, SEXP nodes, SEXP nid, SEXP sublist, SEXP parents);
SEXP arcs2walist(SEXP arcs, SEXP nodes, SEXP weights, SEXP nid, SEXP sublist,
    SEXP parents);

/* convert an arc set to an edge list. */
SEXP arcs2alist(SEXP arcs, SEXP nodes, SEXP weights, SEXP nid, SEXP sublist,
    SEXP parents) {

  if (weights == R_NilValue)
    return arcs2ualist(arcs, nodes, nid, sublist, parents);
  else
    return arcs2walist(arcs, nodes, weights, nid, sublist, parents);

}/*ARCS2ALIST*/

/* convert an arc set to a weighted edge list. */
SEXP arcs2walist(SEXP arcs, SEXP nodes, SEXP weights, SEXP nid, SEXP sublist,
    SEXP parents) {

int i = 0, j = 0, k = 0, nnodes = length(nodes), narcs = length(arcs)/2;
int *e = NULL, *coords = NULL, *adjacent = NULL;
int num_id = isTRUE(nid), sub = isTRUE(sublist), up = isTRUE(parents);
double *w = REAL(weights), *ew = NULL;
SEXP try, alist, edges, edge_weights, temp, temp_name = R_NilValue;

  /* allocate the return value. */
  PROTECT(alist = allocVector(VECSXP, nnodes));
  /* set the node names. */
  setAttrib(alist, R_NamesSymbol, nodes);

  /* allocate and initialize the subset name. */
  if (sub > 0)
    PROTECT(temp_name = mkStringVec(2, "edges", "weight"));

  /* allocate the scratch space to keep track adjacent nodes. */
  adjacent = Calloc1D(nnodes, sizeof(int));

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  for (i = 0; i < narcs; i++)
    adjacent[coords[i + up * narcs] - 1]++;

  for (i = 0; i < nnodes; i++) {

    /* allocate and set up the edge array. */
    if (num_id > 0) {

      PROTECT(edges = allocVector(INTSXP, adjacent[i]));
      e = INTEGER(edges);

    }/*THEN*/
    else {

      PROTECT(edges = allocVector(STRSXP, adjacent[i]));

    }/*ELSE*/

    /* allocate and set up the weights array. */
    PROTECT(edge_weights = allocVector(REALSXP, adjacent[i]));
    ew = REAL(edge_weights);

    /* copy the coordinates or the labels of adjacent nodes. */
    for (j = 0, k = 0; j < narcs; j++) {

      if (coords[j + up * narcs] != i + 1)
        continue;

      /* copy the weight as well. */
      ew[k] = w[j];

      if (num_id > 0)
        e[k++] = coords[(1 - up) * narcs + j];
      else
        SET_STRING_ELT(edges, k++, STRING_ELT(arcs, (1 - up) * narcs + j));

      if (k == adjacent[i])
        break;

    }/*FOR*/

    if (sub > 0) {

      /* allocate and set up the "edge" sublist for graphNEL. */
      PROTECT(temp = allocVector(VECSXP, 2));
      setAttrib(temp, R_NamesSymbol, temp_name);
      SET_VECTOR_ELT(temp, 0, edges);
      SET_VECTOR_ELT(temp, 1, edge_weights);
      SET_VECTOR_ELT(alist, i, temp);
      UNPROTECT(1);

    }/*THEN*/
    else {

      /* save weights with edges as names. */
      setAttrib(edge_weights, R_NamesSymbol, edges);
      SET_VECTOR_ELT(alist, i, edge_weights);

    }/*ELSE*/

    UNPROTECT(2);

  }/*FOR*/

  Free1D(adjacent);

  if (sub > 0)
    UNPROTECT(3);
  else
    UNPROTECT(2);

  return alist;

}/*ARCS2WALIST*/

/* convert an arc set to an unweighted edge list. */
SEXP arcs2ualist(SEXP arcs, SEXP nodes, SEXP nid, SEXP sublist, SEXP parents) {

int i = 0, j = 0, k = 0, nnodes = length(nodes), narcs = length(arcs)/2;
int *e = NULL, *coords = NULL, *adjacent = NULL, up = isTRUE(parents);
int num_id = isTRUE(nid), sub = isTRUE(sublist);
SEXP try, alist, edges, temp, temp_name = R_NilValue;

  /* allocate the return value. */
  PROTECT(alist = allocVector(VECSXP, nnodes));
  /* set the node names. */
  setAttrib(alist, R_NamesSymbol, nodes);

  /* allocate and initialize the subset name. */
  if (sub > 0)
    PROTECT(temp_name = mkString("edges"));

  /* allocate the scratch space to keep track of adjacent nodes. */
  adjacent = Calloc1D(nnodes, sizeof(int));

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  for (i = 0; i < narcs; i++)
    adjacent[coords[i + up * narcs] - 1]++;

  for (i = 0; i < nnodes; i++) {

    /* allocate and set up the edge array. */
    if (num_id > 0) {

      PROTECT(edges = allocVector(INTSXP, adjacent[i]));
      e = INTEGER(edges);

    }/*THEN*/
    else {

      PROTECT(edges = allocVector(STRSXP, adjacent[i]));

    }/*ELSE*/

    /* copy the coordinates or the labels of adjacent nodes. */
    for (j = 0, k = 0; j < narcs; j++) {

      if (coords[j + up * narcs] != i + 1)
        continue;

      if (num_id > 0)
        e[k++] = coords[(1 - up) * narcs + j];
      else
        SET_STRING_ELT(edges, k++, STRING_ELT(arcs, (1 - up) * narcs + j));

      if (k == adjacent[i])
        break;

    }/*FOR*/

    if (sub > 0) {

      /* allocate and set up the "edge" sublist for graphNEL. */
      PROTECT(temp = allocVector(VECSXP, 1));
      setAttrib(temp, R_NamesSymbol, temp_name);
      SET_VECTOR_ELT(temp, 0, edges);
      SET_VECTOR_ELT(alist, i, temp);
      UNPROTECT(1);

    }/*THEN*/
    else {

      SET_VECTOR_ELT(alist, i, edges);

    }/*ELSE*/

    UNPROTECT(1);

  }/*FOR*/

  Free1D(adjacent);

  if (sub > 0)
    UNPROTECT(3);
  else
    UNPROTECT(2);

  return alist;

}/*ARCS2UALIST*/

/* convert an edge list to an arc set. */
SEXP alist2arcs(SEXP alist, SEXP parents) {

int i = 0, j = 0, k = 0, n = length(alist), narcs = 0, up = isTRUE(parents);
SEXP from, adjacent, nodes, arcs;

  /* count how many arcs are present in the graph. */
  for (i = 0; i < n; i++)
    narcs += length(VECTOR_ELT(alist, i));

  /* allocate and initialize the return value. */
  PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
  setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));

  nodes = getAttrib(alist, R_NamesSymbol);

  for (i = 0; i < n; i++) {

     /* cache the parent node and the vector of its adjacent. */
     from = STRING_ELT(nodes, i);
     adjacent = VECTOR_ELT(alist, i);

     /* fill the return value. */
     for (j = 0; j < length(adjacent); j++) {

       SET_STRING_ELT(arcs, k + narcs * up, from);
       SET_STRING_ELT(arcs, k + narcs * (1 - up), STRING_ELT(adjacent, j));
       k++;

     }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);

  return arcs;

}/*ALIST2ARCS*/
