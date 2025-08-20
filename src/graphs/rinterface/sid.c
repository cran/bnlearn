#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../include/graph.h"
#include "../../include/globals.h"
#include "../../minimal/common.h"
#include "../../minimal/strings.h"
#include "../../math/linear.algebra.h"
#include "../graphs.h"

/* build a reachability matrix that takes non-causal paths into account. */
SEXP reachability_matrix(SEXP dag, SEXP path_matrix, SEXP target_node,
    SEXP adjustment_set, SEXP debug) {

int nnodes = 0, target = 0;
int *initial_rmat = NULL, *final_rmat = NULL;
char **labels = NULL;
bool *aset = NULL;
SEXP nodes, arcs, amat, names, reachability, temp;

  /* dereference the node labels and the arcs. */
  arcs = getListElement(dag, "arcs");
  nodes = getListElement(dag, "nodes");
  nodes = getAttrib(nodes, R_NamesSymbol);
  nnodes = length(nodes);
  /* for the node labels, duplicate them for use in the reachability matrix. */
  labels = Calloc1D(2 * nnodes, sizeof(char *));
  for (int i = 0; i < nnodes; i++)
    labels[i] = labels[i + nnodes] = (char *) CHAR(STRING_ELT(nodes, i));

  /* contruct the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));

  /* dereference the target node. */
  PROTECT(temp = match(nodes, target_node, 0));
  target = INT(temp) - 1;
  UNPROTECT(1);

  /* explode the adjustment set into a logical vector to reduce the amount of
   * nested loops later. */
  PROTECT(temp = match(nodes, adjustment_set, 0));
  aset = Calloc1D(nnodes, sizeof(bool));
  for (int i = 0; i < length(adjustment_set); i++)
    aset[INTEGER(temp)[i] - 1] = TRUE;
  UNPROTECT(1);

  /* initialize the reachability matrix. */
  initial_rmat = Calloc1D(4 * nnodes * nnodes, sizeof(int));

  initialize_reachability_matrix(INTEGER(amat), INTEGER(path_matrix), nnodes,
      target, aset, initial_rmat, labels, isTRUE(debug));

  /* allocate the final reachability matrix with all the attributes. */
  PROTECT(reachability = allocMatrix(LGLSXP, 2 * nnodes, 2 * nnodes));
  final_rmat = LOGICAL(reachability);
  memset(final_rmat, '\0', 4 * nnodes * nnodes * sizeof(int));

  PROTECT(names = allocVector(STRSXP, 2 * nnodes));
  for (int i = 0; i < 2 * nnodes; i++)
    SET_STRING_ELT(names, i, STRING_ELT(nodes, i % nnodes));
  setDimNames(reachability, names, names);

  /* traverse the paths in the reachability matrix to complete it. */
  complete_reachability_matrix(initial_rmat, final_rmat, 2 * nnodes, labels,
    isTRUE(debug));

  Free1D(aset);
  Free1D(labels);
  Free1D(initial_rmat);

  UNPROTECT(3);

  return reachability;

}/*REACHABILITY_MATRIX*/
