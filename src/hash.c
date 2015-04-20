#include "include/rcore.h"
#include "include/matrix.h"

SEXP arc_hash_matrix(SEXP arcs, SEXP nodes, int uptri, int sort);
SEXP arc_hash_dataframe(SEXP arcs, SEXP nodes, int uptri, int sort);

/* compute the hash of an arc set using CMC()/UPTRI3() coordinates. */
SEXP arc_hash(SEXP arcs, SEXP nodes, int uptri, int sort) {

  switch(TYPEOF(arcs)) {

    case STRSXP:
      return arc_hash_matrix(arcs, nodes, uptri, sort);
      break;

    case VECSXP:
      return arc_hash_dataframe(arcs, nodes, uptri, sort);
      break;

  }/*SWITCH*/

  return R_NilValue;

}/*ARC_HASH*/

/* arc hash backend. */
void c_arc_hash(int narcs, int nnodes, int *from, int *to, int *uptri,
  int *cmc, int sort) {

int i = 0;

  if (uptri) {

    for (i = 0; i < narcs; i++)
      uptri[i] = UPTRI3(from[i], to[i], nnodes);

    if (sort)
      R_isort(uptri, narcs);

  }/*THEN*/

  if (cmc) {

    for (i = 0; i < narcs; i++)
      cmc[i] = CMC(from[i], to[i], nnodes);

    if (sort)
      R_isort(cmc, narcs);

  }/*ELSE*/

}/*C_ARC_HASH*/

/* compute the hash of an arc set from an arc matrix. */
SEXP arc_hash_matrix(SEXP arcs, SEXP nodes, int uptri, int sort) {

int narcs = length(arcs)/2, nnodes = length(nodes);
int *temp = NULL, *hash_el = NULL;
SEXP matched, hash;

  /* match the lables of the nodes. */
  PROTECT(matched = match(nodes, arcs, 0));
  temp = INTEGER(matched);

  /* allocate the integer array to hold the hash. */
  PROTECT(hash = allocVector(INTSXP, narcs));
  hash_el = INTEGER(hash);

  if (uptri)
    c_arc_hash(narcs, nnodes, temp, temp + narcs, hash_el, NULL, sort);
  else
    c_arc_hash(narcs, nnodes, temp, temp + narcs, NULL, hash_el, sort);

  UNPROTECT(2);

  return hash;

}/*ARC_HASH_MATRIX*/

/* compute the hash of an arc set from a data frame. */
SEXP arc_hash_dataframe(SEXP arcs, SEXP nodes, int uptri, int sort) {

int narcs = length(VECTOR_ELT(arcs, 0)), nnodes = length(nodes);
int *temp1 = NULL, *temp2 = NULL, *hash_el = NULL;
SEXP matched1, matched2, hash;

  /* match the labels in the first column. */
  PROTECT(matched1 = match(nodes, VECTOR_ELT(arcs, 0), 0));
  temp1 = INTEGER(matched1);
  /* match the labels in the second column. */
  PROTECT(matched2 = match(nodes, VECTOR_ELT(arcs, 1), 0));
  temp2 = INTEGER(matched2);

  /* allocate the integer array to hold the hash. */
  PROTECT(hash = allocVector(INTSXP, narcs));
  hash_el = INTEGER(hash);

  if (uptri)
    c_arc_hash(narcs, nnodes, temp1, temp2, hash_el, NULL, sort);
  else
    c_arc_hash(narcs, nnodes, temp1, temp2, NULL, hash_el, sort);

  UNPROTECT(3);

  return hash;

}/*ARC_HASH_DATAFRAME*/

/* compute the hash of an adjacency matrix using CMC() coordinates. */
SEXP c_amat_hash(int *amat, int nnodes) {

int i = 0, k = 0, narcs = 0;
int *coords = NULL;
SEXP hash;

  /* count the arcs in the network.  */
  for (i = 0; i < nnodes * nnodes; i++)
    if (amat[i] > 0)
      narcs++;

  /* allocate the hash (a vector containing the coordinates of the arcs in
   * the flattened adjacency matrix, which uniquely identifies a directed
   * graph). */
  PROTECT(hash = allocVector(INTSXP, narcs));
  coords = INTEGER(hash);

  /* the (flattened) coordinates into the hash vector. */
  for (i = 0, k = 0; i < nnodes * nnodes; i++)
    if (amat[i] > 0)
      coords[k++] = i;

  UNPROTECT(1);

  return hash;

}/*C_AMAT_HASH*/

