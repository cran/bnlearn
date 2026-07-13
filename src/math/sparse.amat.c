#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/sort.h"
#include "../minimal/common.h"
#include "sparse.amat.h"

/* build a CSR adjacency matrix from an arc set. */
sparse_amat new_sparse_amat(int *from, int *to, int narcs, int nnodes) {

sparse_amat sam = { 0 };
int *pos = NULL;

  sam.dim = nnodes;
  sam.narcs = narcs;

  /* count the children of each node to lay out the per-node offsets. */
  sam.rowptr = Calloc1D(nnodes + 1, sizeof(int));
  for (int k = 0; k < narcs; k++)
    sam.rowptr[from[k] + 1]++;
  for (int i = 0; i < nnodes; i++)
    sam.rowptr[i + 1] += sam.rowptr[i];

  /* place each child into its parent's slice of the children array. */
  sam.colidx = Calloc1D(narcs, sizeof(int));
  pos = Calloc1D(nnodes, sizeof(int));
  for (int i = 0; i < nnodes; i++)
    pos[i] = sam.rowptr[i];
  for (int k = 0; k < narcs; k++)
    sam.colidx[pos[from[k]]++] = to[k];
  Free1D(pos);

  /* sort each node's children by index, to match the order in which a dense
   * adjacency matrix would have been scanned. */
  for (int i = 0; i < nnodes; i++)
    i_sort(sam.colidx + sam.rowptr[i], NULL, sam.rowptr[i + 1] - sam.rowptr[i]);

  return sam;

}/*NEW_SPARSE_AMAT*/

/* build a CSR adjacency matrix from a "bn" object: match its arcs against its own
 * node labels and call new_sparse_amat(), optionally transposing (so that the
 * children matrix becomes the parents matrix), then attach the node labels. */
sparse_amat sparse_amat_from_SEXP(SEXP bn, bool transpose) {

sparse_amat sam = { 0 };
int nnodes = 0, narcs = 0;
int *from = NULL, *to = NULL;
SEXP arcs, nodes, arc_ids;

  arcs = getListElement(bn, "arcs");
  PROTECT(nodes = getAttrib(getListElement(bn, "nodes"), R_NamesSymbol));
  nnodes = length(nodes);
  narcs = length(arcs) / 2;

  /* map the arc endpoints to node indexes. */
  PROTECT(arc_ids = match(nodes, arcs, 0));
  from = Calloc1D(narcs, sizeof(int));
  to = Calloc1D(narcs, sizeof(int));
  for (int k = 0; k < narcs; k++) {

    from[k] = INTEGER(arc_ids)[k] - 1;
    to[k] = INTEGER(arc_ids)[k + narcs] - 1;

  }/*FOR*/

  /* transposing swaps the arc endpoints, turning children into parents. */
  if (transpose)
    sam = new_sparse_amat(to, from, narcs, nnodes);
  else
    sam = new_sparse_amat(from, to, narcs, nnodes);

  Free1D(from);
  Free1D(to);

  /* attach the row/column (node) labels. */
  sam.labels = Calloc1D(nnodes, sizeof(char *));
  for (int k = 0; k < nnodes; k++)
    sam.labels[k] = (char *) CHAR(STRING_ELT(nodes, k));

  UNPROTECT(2);

  return sam;

}/*SPARSE_AMAT_FROM_SEXP*/

void FreeSparseAMAT(sparse_amat sam) {

  Free1D(sam.rowptr);
  Free1D(sam.colidx);
  Free1D(sam.labels);

}/*FREESPARSEAMAT*/
