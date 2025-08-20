#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../math/linear.algebra.h"

/* build a path matrix, telling whether each node is reachable from another. */
void dag_path_matrix(int *amat, int *pathmat, int nnodes, char **labels,
    int *root_ids, int nroots, bool debugging) {

int depth = 0, cur = 0, root = 0;
int *counters = NULL, *path = NULL;
bool *visited = NULL, moved = FALSE;

  /* status vector for the nodes. */
  visited = Calloc1D(nnodes, sizeof(bool));
  /* initialize the position counters for the rows of the adjacency matrix. */
  counters = Calloc1D(nnodes, sizeof(int));
  /* allocate the path array. */
  path = Calloc1D(nnodes, sizeof(int));

  /* for each node... */
  for (int node = 0; node < nroots; node++) {

    root = root_ids[node] - 1;

    if (debugging)
      Rprintf("* checking node %s.\n", labels[root]);

    /* ... reset the position counters... */
    memset(counters, '\0', nnodes * sizeof(int));

    /* ... initialize the path array for that node. */
    memset(path, '\0', nnodes * sizeof(int));
    path[0] = root;
    depth = 0;

    cur = root;

    /* until we have explored all possible paths with a depth-first search and
     * back-tracked all the way to where we started... */
    do {

      moved = FALSE;

      /* ... find the next child of the last node in the path... */
      for (int i = counters[cur]; i < nnodes; i++) {

        counters[cur]++;

        if (amat[CMC(cur, i, nnodes)] > 0) {

          /* ... and add it to the path. */
          path[++depth] = i;
          cur = i;

          /* the last node is reachable from all earlier nodes in the path. */
          for (int k = 0; k < depth; k++)
            pathmat[CMC(path[k], path[depth], nnodes)] = TRUE;

          if (debugging) {

            Rprintf("  > found node %s.\n", labels[i]);
            Rprintf("  > path: ");
            for (int k = 0; k < depth; k++)
              Rprintf("%s -> ", labels[path[k]]);
            Rprintf("%s\n", labels[path[depth]]);

          }/*THEN*/

          moved = TRUE;
          break;

        }/*THEN*/

      }/*FOR*/

      /* check whether to update the path matrix or backtrack. */
      if (!visited[path[depth]]) {

        visited[path[depth]] = TRUE;

      }/*THEN*/
      else {

        if (moved) {

          if (debugging)
            Rprintf("  @ already visited node %s.\n", labels[path[depth]]);

          /* we already know which nodes are reachable from here; update the
           * path matrix ... */
          for (int j = 0; j < nnodes; j++) {

            if (!pathmat[CMC(path[depth], j, nnodes)])
              continue;

            for (int k = 0; k < depth; k++) {

              pathmat[CMC(path[k], j, nnodes)] = TRUE;

              if (debugging)
                Rprintf("  & node %s is reachable from %s, so it is also reachable from node %s.\n",
                    labels[j], labels[path[depth]], labels[path[k]]);

            }/*FOR*/

          }/*FOR*/

        }/*THEN*/

        /* ... and backtrack. */
        if (depth > 0) {

          cur = path[--depth];
          path[depth + 1] = 0;

        }/*THEN*/

      }/*ELSE*/

    } while (!((depth == 0) && (counters[root] == nnodes)));

  }/*FOR*/

  Free1D(counters);
  Free1D(path);
  Free1D(visited);

  /* complete the path matrix by setting the diagonal. */
  for (int i = 0; i < nnodes; i++)
    pathmat[CMC(i, i, nnodes)] = TRUE;

}/*DAG_PATH_MATRIX*/

