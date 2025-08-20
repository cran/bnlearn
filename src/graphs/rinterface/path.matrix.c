#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../include/graph.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/strings.h"
#include "../../math/linear.algebra.h"
#include "../graphs.h"

/* build a path matrix, telling whether each node is reachable from another. */
SEXP path_matrix(SEXP x, SEXP debug) {

int nnodes = 0;
char **labels = NULL;
SEXP arcs, nodes, amat, roots, root_ids, pathmat;

  /* extract the relevant information from the network. */
  arcs = getListElement(x, "arcs");
  nodes = getListElement(x, "nodes");
  nodes = getAttrib(nodes, R_NamesSymbol);
  nnodes = length(nodes);
  /* dereference the node labels. */
  labels = Calloc1D(nnodes, sizeof(char *));
  for (int i = 0; i < nnodes; i++)
    labels[i] = (char *) CHAR(STRING_ELT(nodes, i));

  /* identify the root nodes. */
  PROTECT(roots = root_nodes(x, FALSESEXP));
  PROTECT(root_ids = match(nodes, roots, 0));

  /* contruct the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));

  /* allocate the path matrix. */
  PROTECT(pathmat = allocMatrix(LGLSXP, nnodes, nnodes));
  setAttrib(pathmat, R_DimNamesSymbol, getAttrib(amat, R_DimNamesSymbol));
  memset(INTEGER(pathmat), '\0', nnodes * nnodes * sizeof(int));

  dag_path_matrix(INTEGER(amat), LOGICAL(pathmat), nnodes, labels,
    INTEGER(root_ids), length(root_ids), isTRUE(debug));

  UNPROTECT(4);

  Free1D(labels);

  return pathmat;

}/*PATH_MATRIX*/

