#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../math/linear.algebra.h"

/* identify the connected components in an undirected graph. */
int ug_connected_components(int *amat, char **labels, int nnodes, int **buffer,
    int *buflen, bool debugging) {

int *elements = NULL, nelements = 0, ncomponents = 0;
bool changed = TRUE, *visited = NULL;

  /* status vector for the nodes. */
  visited = Calloc1D(nnodes, sizeof(bool));
  /* depth of the nodes in the breadth-first search. */
  elements = Calloc1D(nnodes, sizeof(int));

  for (int i = 0; i < nnodes; i++) {

    /* for each node we have not visited already... */
    if (visited[i])
      continue;
    else if (debugging)
      Rprintf("* checking node %s.\n", labels[i]);

    /* ... mark it as visited. */
    visited[i] = TRUE;

    /* reset the connected component... */
    memset(elements, '\0', nnodes * sizeof(int));
    nelements = 1;
    elements[i] = 1;
    changed = FALSE;

    /* in a breadth-first search... */
    for (int depth = 2; depth < nnodes; depth++) {

      /* ... check the remaining nodes... */
      for (int j = 0; j < nnodes; j++) {

        /* .. if we have visited them in the last iteration... */
        if (elements[j] != depth - 1)
          continue;

        /* ... mark their neighbours. */
        for (int k = 0; k < nnodes; k++) {

          if ((amat[CMC(j, k, nnodes)] != 0) && (elements[k] == 0)) {

            elements[k] = depth;
            visited[k] = TRUE;
            changed = TRUE;
            nelements++;

          }/*THEN*/

        }/*FOR*/

      }/*FOR*/

      /* if we found no new nodes, we are done with this component. */
      if (changed)
        changed = FALSE;
      else
        break;

    }/*FOR*/

    /* save the elements of the connected component. */
    buffer[ncomponents] = Calloc1D(nelements, sizeof(int));
    buflen[ncomponents] = nelements;
    for (int k = 0, l = 0; k < nnodes; k++)
      if (elements[k] > 0)
        buffer[ncomponents][l++] = k;

    if (debugging) {

      Rprintf("  @ connected component %d, %d nodes:\n",
        ncomponents, nelements);
      for (int k = 0; k < buflen[ncomponents]; k++)
        Rprintf("%s ", labels[buffer[ncomponents][k]]);
      Rprintf("\n");

    }/*THEN*/

    /* move to the next component. */
    ncomponents++;

  }/*FOR*/

  Free1D(visited);
  Free1D(elements);

  return ncomponents;

}/*UG_CONNECTED_COMPONENTS*/

