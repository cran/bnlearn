#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/sort.h"
#include "fitted.h"

/* sort the nodes of a BN in topological order. */
void topological_sort_bn(fitted_bn fitted, int *poset) {

int nnodes = fitted.nnodes, *depth = NULL;
bool changed = TRUE;

  depth = Calloc1D(nnodes, sizeof(int));

  /* root nodes (no parents) are at depth 1, all others unassigned (0). */
  for (int i = 0; i < nnodes; i++)
    depth[i] = (fitted.ldists[i].nparents == 0) ? 1 : 0;

  /* relax the depths until they stabilize. */
  while (changed) {

    changed = FALSE;

    for (int i = 0; i < nnodes; i++) {

      int np = fitted.ldists[i].nparents, maxdepth = 0, ready = TRUE;

      if (np == 0)
        continue;

      for (int k = 0; k < np; k++) {

        int pd = depth[fitted.ldists[i].parents[k]];

        if (pd == 0) {

          ready = FALSE;
          break;

        }/*THEN*/

        if (pd > maxdepth)
          maxdepth = pd;

      }/*FOR*/

      if (ready && (depth[i] != maxdepth + 1)) {

        depth[i] = maxdepth + 1;
        changed = TRUE;

      }/*THEN*/

    }/*FOR*/

  }/*WHILE*/

  /* sort the node indexes by their depth. */
  for (int i = 0; i < nnodes; i++)
    poset[i] = i;
  i_sort(depth, poset, nnodes);

  Free1D(depth);

}/*TOPOLOGICAL_SORT_BN*/
