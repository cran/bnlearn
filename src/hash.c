#include "common.h"

/* compute the hash of an arc set using CMC() coordinates. */
SEXP arc_hash(SEXP arcs, SEXP nodes) {

int i = 0, narcs = LENGTH(arcs)/2;
int *temp = NULL, *hash_el = NULL;
SEXP matched, hash;

  /* match the lables of the nodes. */
  PROTECT(matched = match(nodes, arcs, 0));
  temp = INTEGER(matched);

  /* allocate the integer array to hold the hash. */
  PROTECT(hash = allocVector(INTSXP, narcs));
  hash_el = INTEGER(hash);

  for (i = 0; i < narcs; i++)
    hash_el[i] = CMC(temp[i], temp[i + narcs], narcs);

  UNPROTECT(2);

  return hash;

}/*ARC_HASH*/

/* compute the hash of an adjacency matrix using CMC() coordinates. */
SEXP c_amat_hash(int *amat, int *nnodes) {

int i = 0, k = 0, narcs = 0;
int *coords = NULL;
SEXP hash;

  /* count the arcs in the network.  */
  for (i = 0; i < (*nnodes) * (*nnodes); i++)
    if (amat[i] > 0)
      narcs++;

  /* allocate the hash (a vector containing the coordinates of the arcs in
   * the flattened adjacency matrix, which uniquely identifies a directed
   * graph). */
  PROTECT(hash = allocVector(INTSXP, narcs));
  coords = INTEGER(hash);

  /* the (flattened) coordinates into the hash vector. */
  for (i = 0, k = 0; i < (*nnodes) * (*nnodes); i++)
    if (amat[i] > 0)
      coords[k++] = i;

  UNPROTECT(1);

  return hash;

}/*C_AMAT_HASH*/

