#include "common.h"

SEXP smart_network_averaging(SEXP arcs, SEXP nodes, SEXP weights) {

int k = 0, from = 0, to = 0, nrows = LENGTH(arcs) / 2, dims = LENGTH(nodes);
int *a = NULL, *coords = NULL, *poset = NULL;
double *w = NULL;
SEXP weights2, amat, try, acyclic;

  /* allocate and initialize the adjacency matrix. */
  PROTECT(amat = allocMatrix(INTSXP, dims, dims));
  a = INTEGER(amat);
  memset(a, '\0', sizeof(int) * dims * dims);

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  /* duplicate the weights to preserve the oroginal ones. */
  PROTECT(weights2 = duplicate(weights));
  w = REAL(weights2);

  /* sort the strength coefficients. */
  poset = alloc1dcont(nrows);
  for (k = 0; k < nrows; k++)
    poset[k] = k;
  R_qsort_I(w, poset, 1, nrows);

  /* iterate over the arcs in reverse order wrt their strength coefficients. */
  for (k = 0; k < nrows; k++) {

    from = coords[poset[k]] - 1;
    to = coords[poset[k] + nrows] - 1;

    /* add an arc only if it does not introduce cycles. */
    if (!c_has_path(to, from, a, dims, nodes, FALSE, TRUE, FALSE))
      a[CMC(from, to, dims)] = 1;
    else
      warning("arc %s -> %s would introduce cycles in the graph, ignoring.",
        NODE(from), NODE(to));

  }/*FOR*/

  /* convert the adjacency matrix back to an arc set and return it. */
  acyclic = amat2arcs(amat, nodes);

  UNPROTECT(3);

  return acyclic;

}/*SMART_NETWORK_AVERAGING*/

