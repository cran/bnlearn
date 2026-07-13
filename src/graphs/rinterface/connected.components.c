#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../include/graph.h"
#include "../../minimal/common.h"
#include "../graphs.h"

/* identify the connected components in an undirected graph. */
SEXP connected_components(SEXP x, SEXP debug) {

int nnodes = 0, ncomponents = 0;
int **buffer = NULL, *buflen = NULL;
char **labels = NULL;
SEXP arcs, nodes, amat, temp, components;

  /* extract the relevant information from the network. */
  arcs = getListElement(x, "arcs");
  nodes = getListElement(x, "nodes");
  PROTECT(nodes = getAttrib(nodes, R_NamesSymbol));
  nnodes = length(nodes);

  /* contruct the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));

  /* buffer to store the connected components. */
  buffer = Calloc1D(nnodes, sizeof(int **));
  buflen = Calloc1D(nnodes, sizeof(int));
  /* dereference the node labels. */
  labels = Calloc1D(nnodes, sizeof(char *));
  for (int i = 0; i < nnodes; i++)
    labels[i] = (char *) CHAR(STRING_ELT(nodes, i));

  /* identify the connected components. */
  ncomponents = ug_connected_components(INTEGER(amat), labels, nnodes, buffer,
                  buflen, isTRUE(debug));

  /* return the components and the respective sets of node labels. */
  PROTECT(components = allocVector(VECSXP, ncomponents));

  for (int i = 0; i < ncomponents; i++) {

    PROTECT(temp = allocVector(STRSXP, buflen[i]));
    for (int j = 0; j < buflen[i]; j++)
      SET_STRING_ELT(temp, j, STRING_ELT(nodes, buffer[i][j]));
    SET_VECTOR_ELT(components, i, temp);

    UNPROTECT(1);

  }/*FOR*/

  UNPROTECT(3);

  Free2D(buffer, nnodes);
  Free1D(buflen);
  Free1D(labels);

  return components;

}/*CONNECTED_COMPONENTS*/
