#include "../include/rcore.h"
#include "../include/graph.h"
#include "../include/globals.h"
#include "../minimal/common.h"
#include "../core/sort.h"

/* find out the partial ordering of the nodes of a DAG. */
SEXP topological_ordering(SEXP bn, SEXP root_nodes, SEXP reverse, SEXP debug) {

int *depth = NULL, *matched = NULL;
int d = 0, i = 0, j = 0, changed = 0, nnodes = 0;
char *direction = NULL;
bool debugging = isTRUE(debug);
SEXP nodes_data, nodes, try, children, ordering;

  if (isTRUE(reverse))
    direction = "parents";
  else
    direction = "children";

  /* get to the nodes' data in both 'bn' and 'bn.fit' objects. */
  nodes_data = getListElement(bn, "nodes");

  if (isNull(nodes_data))
    nodes_data = bn;

  /* get and count the node labels. */
  PROTECT(nodes = getAttrib(nodes_data, R_NamesSymbol));
  nnodes = length(nodes);

  /* allocate a status vector to track the ordering of the nodes. */
  PROTECT(ordering = allocVector(INTSXP, nnodes));
  depth = INTEGER(ordering);
  memset(depth, '\0', nnodes * sizeof(int));

  if (debugging)
    Rprintf("* currently at depth 1 (starting BFS).\n");

  /* set the root nodes as the starting point of the BFS. */
  PROTECT(try = match(nodes, root_nodes, 0));
  matched = INTEGER(try);

  for (i = 0; i < length(try); i++) {

    if (debugging)
      Rprintf("  > got node %s.\n", NODE(matched[i] - 1));

    depth[matched[i] - 1] = 1;

  }/*FOR*/

  UNPROTECT(1);

  /* now let's go down from the roots to the leafs, one layer at a time. */
  for (d = 1; d <= nnodes; d++) {

    if (debugging)
      Rprintf("* currently at depth %d.\n", d + 1);

    /* reset the changed flag. */
    changed = 0;

    for (i = 0; i < nnodes; i++) {

      /* this node has already been visisted, skip. */
      if (depth[i] < d)
        continue;

      children = getListElement(VECTOR_ELT(nodes_data, i), direction);

      /* this node is a leaf, nothing to do, move along. */
      if (length(children) == 0)
        continue;

      /* set the changed flag. */
      changed = 1;

      PROTECT(try = match(nodes, children, 0));
      matched = INTEGER(try);

      /* set the correct depth to the children of this node. */
      for (j = 0; j < length(try); j++) {

        if (debugging)
          Rprintf("  > got node %s from %s.\n",
            NODE(matched[j] - 1), NODE(i));

        depth[matched[j] - 1] = d + 1;

      }/*FOR*/

      UNPROTECT(1);

    }/*FOR*/

    /* all nodes have been visited, break. */
    if (!changed) break;

  }/*FOR*/

  if (debugging > 0)
    Rprintf("* all nodes have been visited.\n");

  /* add the node labels to the return value. */
  setAttrib(ordering, R_NamesSymbol, nodes);

  UNPROTECT(2);

  return ordering;

}/*TOPOLOGICAL_ORDERING*/

/* sort the nodes in topological order. */
void topological_sort(SEXP fitted, int *poset, int nnodes) {

int i = 0;
SEXP roots, node_depth;

  PROTECT(roots = root_nodes(fitted, FALSESEXP));
  PROTECT(node_depth = topological_ordering(fitted, roots, FALSESEXP, FALSESEXP));
  for (i = 0; i < nnodes; i++)
    poset[i] = i;
  i_sort(INTEGER(node_depth), poset, nnodes);

  UNPROTECT(2);

}/*TOPOLOGICAL_SORT*/
