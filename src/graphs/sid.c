#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../math/linear.algebra.h"

/* small helper function to avoid duplicating code everywhere. */
static bool already_scheduled(int node, char *label, int *set, int size,
    bool debugging) {

bool found = FALSE;

  for (int already = 0; already < size; already++)
    if (set[already] == node) {

      found = TRUE;
      break;

    }/*THEN*/

  if (debugging) {

    if (found)
      Rprintf("    @ node %s already scheduled for checking.\n", label);
    else
      Rprintf("    @ schedule node %s for checking.\n", label);

  }/*THEN*/

  return found;

}/*ALREADY_SCHEDULED*/

/* initialize a reachability matrix that takes non-causal paths into account. */
void initialize_reachability_matrix(int *amat, int *pmat, int nnodes, int target,
    bool *adjustment_set, int *rmat, char **labels, bool debugging) {

int *to_check = NULL, last_node = 0;
bool *checked = NULL, *adjustment_ancestors = NULL;

  /* allocate the status vectors, assuming we only check each node once. */
  to_check = Calloc1D(nnodes, sizeof(int));
  checked = Calloc1D(nnodes, sizeof(bool));

  /* the initial set of nodes to check: the parents and children. */
  for (int i = 0; i < nnodes; i++)
    if ((amat[CMC(i, target, nnodes)] > 0) || (amat[CMC(target, i, nnodes)] > 0))
      to_check[last_node++] = i;

  /* find the ancestors of the adjustment set, and mark the adjustment set as
   * well for unknown reasons. */
  adjustment_ancestors = Calloc1D(nnodes, sizeof(bool));

  for (int i = 0; i < nnodes; i++) {

    if (!adjustment_set[i])
      continue;

    adjustment_ancestors[i] = TRUE;

    for (int j = 0; j < nnodes; j++)
      if (pmat[CMC(j, i, nnodes)])
        adjustment_ancestors[j] = TRUE;

  }/*FOR*/

  /* nuke all parents and children of the target node. */
  for (int i = 0; i < nnodes; i++)
    amat[CMC(i, target, nnodes)] = amat[CMC(target, i, nnodes)] = 0;

  if (debugging)
    Rprintf("* creating the reachability matrix for node %s.\n", labels[target]);

  /* for each node... */
  for (int cur = 0; cur < nnodes; cur++) {

    /* ... assuming that we still have nodes to check... */
    if (cur >= last_node)
      break;

    if (debugging)
      Rprintf("  > checking node %s\n", labels[to_check[cur]]);

    /* if one of the parents of the current node is reachable and is not part of
     * the adjustment set, the current node is reachable too. */
    for (int pa = 0; pa < nnodes; pa++) {

      if (!((amat[CMC(pa, to_check[cur], nnodes)] > 0) && (!adjustment_set[pa])))
        continue;

      if (debugging)
        Rprintf("    > found parent %s not in the adjustment set.\n", labels[pa]);

      rmat[CMC(pa, to_check[cur], 2 * nnodes)] = TRUE;
      rmat[CMC(nnodes + pa, to_check[cur], 2 * nnodes)] = TRUE;

    }/*FOR*/

    /* if the current node is an ancestor of the target node... */
    if (adjustment_ancestors[to_check[cur]]) {

      /* ... its parents are as well. */
      for (int pa = 0; pa < nnodes; pa++) {

        if (amat[CMC(pa, to_check[cur], nnodes)] == 0)
          continue;

        if (debugging)
          Rprintf("    > found parent %s with ancestors.\n", labels[pa]);

        rmat[CMC(to_check[cur], nnodes + pa, 2 * nnodes)] = TRUE;

        /* add the parents to the nodes that will be checked. */
        if (!already_scheduled(pa, labels[pa], to_check, last_node, debugging))
          to_check[last_node++] = pa;

      }/*FOR*/

    }/*THEN*/

    /* the current node is not in the adjustment set... */
    if (!adjustment_set[to_check[cur]]) {

      /* ... mark the parents as reachable, too. */
      for (int pa = 0; pa < nnodes; pa++) {

        if (amat[CMC(pa, to_check[cur], nnodes)] == 0)
          continue;

        if (debugging)
          Rprintf("    > found parent %s while checking node not in adjustment set.\n",
              labels[pa]);

        rmat[CMC(nnodes + to_check[cur], nnodes + pa, 2 * nnodes)] = TRUE;

        /* add the parents to the nodes that will be checked. */
        if (!already_scheduled(pa, labels[pa], to_check, last_node, debugging))
          to_check[last_node++] = pa;

      }/*FOR*/

    }/*THEN*/

    /* if a child of the current node is reachable and not in the adjustment
     * set, then the current node is too. */
    for (int ch = 0; ch < nnodes; ch++) {

      if (!((amat[CMC(to_check[cur], ch, nnodes)] > 0) && (!adjustment_set[ch])))
        continue;

      if (debugging)
        Rprintf("    > found child %s not in the adjustment set.\n", labels[ch]);

      rmat[CMC(nnodes + ch, nnodes + to_check[cur], 2 * nnodes)] = TRUE;

    }/*FOR*/

    /* if a child of the current node is reachable and and it is an ancestor of
     * the target node, then the current node is reachable. */
    for (int ch = 0; ch < nnodes; ch++) {

      if (!((amat[CMC(to_check[cur], ch, nnodes)] > 0) && (adjustment_ancestors[ch])))
        continue;

      if (debugging)
        Rprintf("    > found child %s among the ancestors.\n", labels[ch]);

      rmat[CMC(ch, nnodes + to_check[cur], 2 * nnodes)] = TRUE;

    }/*FOR*/

    if (!adjustment_set[to_check[cur]]) {

      for (int ch = 0; ch < nnodes; ch++) {

        if (amat[CMC(to_check[cur], ch, nnodes)] == 0)
          continue;

        if (debugging)
          Rprintf("    > found child %s while checking node not in adjustment set.\n",
              labels[ch]);

        rmat[CMC(to_check[cur], ch, 2 * nnodes)] = TRUE;
        rmat[CMC(nnodes + to_check[cur], ch, 2 * nnodes)] = TRUE;

        /* add the children to the nodes that will be checked. */
        if (!already_scheduled(ch, labels[ch], to_check, last_node, debugging))
          to_check[last_node++] = ch;

      }/*FOR*/

    }/*THEN*/

  }/*FOR*/

  Free1D(to_check);
  Free1D(checked);
  Free1D(adjustment_ancestors);

}/*INITIALIZE_REACHABILITY_MATRIX*/

/* complete a reachability matrix by traversing all paths. */
void complete_reachability_matrix(int *initial_rmat, int *final_rmat,
    int nnodes, char **labels, bool debugging) {

int depth = 0, cur = 0, root = 0, ncycles = 0, cycles_capacity = 4 * nnodes;
int *counters = NULL, *path = NULL, **cycles = NULL;
bool *roots = NULL, *completed = NULL, *visited = NULL;
bool moved = FALSE, backtrack = FALSE;

  /* status vectors for the nodes. */
  completed = Calloc1D(nnodes, sizeof(bool));
  visited = Calloc1D(nnodes, sizeof(bool));
  /* initialize the position counters for the rows of the adjacency matrix. */
  counters = Calloc1D(nnodes, sizeof(int));
  /* allocate the path array. */
  path = Calloc1D(nnodes, sizeof(int));
  /* allocate the array to store the cycles. */
  cycles = Calloc1D(cycles_capacity, sizeof(int *));

  /* identify the root nodes. */
  roots = Calloc1D(nnodes, sizeof(int));

  for (int j = 0; j < nnodes; j++) {

    bool found_parent = FALSE;

    for (int i = 0; i < nnodes; i++)
      if (initial_rmat[CMC(i, j, nnodes)]) {

        found_parent = TRUE;
        break;
      }/*THEN*/

    if (found_parent)
      continue;

    roots[j] = TRUE;

  }/*FOR*/

  /* for each node... */
  for (int node = 0; node < nnodes; node++) {

    if (roots[node])
      root = node;
    else
      continue;

    if (debugging)
      Rprintf("  * checking node %s.\n", labels[root]);

    /* ... reset the position counters and the cycle counter... */
    memset(counters, '\0', nnodes * sizeof(int));
    memset(completed, '\0', nnodes * sizeof(bool));
    ncycles = 0;

    /* ... initialize the path array for that node. */
    memset(path, '\0', nnodes * sizeof(int));
    path[0] = root;
    depth = 0;

    cur = root;

    /* until we have explored all possible paths with a depth-first search and
     * back-tracked all the way to where we started... */
    do {

      moved = FALSE;
      backtrack = FALSE;

      /* ... find the next child of the last node in the path... */
      for (int i = counters[cur]; i < nnodes; i++) {

        counters[cur]++;

        if (initial_rmat[CMC(cur, i, nnodes)] > 0) {

          /* ... and add it to the path. */
          path[++depth] = i;
          cur = i;

          /* the last node is reachable from all earlier nodes in the path. */
          for (int k = 0; k < depth; k++)
            final_rmat[CMC(path[k], path[depth], nnodes)] = TRUE;

          /* have already passed through all the arcs going out of this node? */
          if (counters[cur] == nnodes)
            completed[cur] = TRUE;

          if (debugging) {

            Rprintf("    > path: ");
            for (int k = 0; k < depth; k++)
              Rprintf("%s -> ", labels[path[k]]);
            Rprintf("%s", labels[path[depth]]);
            if (visited[cur])
              Rprintf(" (explored all outgoing paths)");
            else if (completed[cur])
              Rprintf(" (explored all outgoing arcs)");
            Rprintf(".\n");

          }/*THEN*/

          moved = TRUE;
          break;

        }/*THEN*/

      }/*FOR*/

      if (!moved) {

        /* if we did not move forward, we can only go back. */
        backtrack = TRUE;

      }/*THEN*/
      else {

        /* cycle detection. */
        for (int j = 0; j < depth; j++)
          if (path[depth] == path[j]) {

            if (debugging)
              Rprintf("    @ node %s already visited in this path, cycle detected.\n",
                  labels[path[depth]]);

            /* allocate space for the cycle + length... */
            if (ncycles >= cycles_capacity) {

              cycles_capacity *= 2;
              cycles = Realloc1D(cycles, cycles_capacity, sizeof(int *));

            }/*THEN*/

            cycles[ncycles] = Calloc1D(depth - j + 1, sizeof(int));
            /* ... copy the node indexes... */
            cycles[ncycles][0] = depth - j;
            memcpy(cycles[ncycles] + 1, path + j, (depth - j) * sizeof(int));
            /* ... and update the counter. */
            ncycles++;

            if (debugging) {

              Rprintf("    [%d]", ncycles - 1);

              for (int j = 0; j < cycles[ncycles - 1][0]; j++)
                Rprintf( " %s", labels[cycles[ncycles - 1][j + 1]]);
              Rprintf("\n");

            }/*FOR*/

            /* backtrack to avoid running in circles. */
            backtrack = TRUE;

          }/*THEN*/

        /* have already passed through all the arcs going out of this node? */
        if (completed[path[depth]] || visited[path[depth]]) {

          /* then we already know all the nodes that are reachable from here,
           * and we can just update the reachability matrix. */
          for (int j = 0; j < nnodes; j++)
            if (final_rmat[CMC(path[depth], j, nnodes)]) {

              for (int k = 0; k < depth; k++)
                final_rmat[CMC(path[k], j, nnodes)] = TRUE;

              if (debugging) {

                Rprintf("    @ %s is then reachable from", labels[j]);
                for (int k = 0; k < depth - 1; k++)
                  Rprintf(" %s,", labels[path[k]]);
                Rprintf(" %s.\n", labels[path[depth - 1]]);

              }/*THEN*/

            }/*THEN*/

          /* all the paths going ouf this node has been completely explored
           * while examining a previous root node, do not explore them again. */
          if (visited[path[depth]])
            backtrack = TRUE;

        }/*THEN*/

      }/*ELSE*/

      if (backtrack && (depth > 0)) {

        cur = path[--depth];
        path[depth + 1] = 0;

      }/*THEN*/

    } while (!((depth == 0) && (counters[root] == nnodes)));

    /* at this point, we have definitely explored all the possible paths going
     * out of the nodes reachable from this particular root node, and we can
     * mark them as visited. */
    for (int j = 0; j < nnodes; j++)
      visited[j] |= completed[j];

    /* adjust for cycles: for each cycle... */
    for (int i = 0; i < ncycles; i++) {

      /* ... for each node in the cycle... */
      for (int j = 0; j < cycles[i][0]; j++)
        /* ... take each other node in the cycle... */
        for (int k = 0; k < cycles[i][0]; k++) {

          /* ... add all its reachable nodes...*/
          for (int l = 0; l < nnodes; l++) {

            if (debugging) {

              if (!final_rmat[CMC(cycles[i][k + 1], l, nnodes)] &&
                   final_rmat[CMC(cycles[i][j + 1], l, nnodes)])
                     Rprintf("  # %s is now reachable from %s.\n",
                       labels[l], labels[cycles[i][k + 1]]);

            }/*THEN*/

            final_rmat[CMC(cycles[i][k + 1], l, nnodes)] |=
              final_rmat[CMC(cycles[i][j + 1], l, nnodes)];

          }/*FOR*/

          /* ... and make those reachable from all the nodes from which we can
           * reach the cycle (reachability is transitive). */
          for (int l = 0; l < nnodes; l++) {

            if (final_rmat[CMC(l, cycles[i][k + 1], nnodes)]) {

              final_rmat[CMC(l, cycles[i][j + 1], nnodes)] |=
                final_rmat[CMC(l, cycles[i][k + 1], nnodes)];

              for (int m = 0; m < nnodes; m++)
                final_rmat[CMC(l, m, nnodes)] |=
                  final_rmat[CMC(cycles[i][k + 1], m, nnodes)];

            }/*THEN*/

          }/*FOR*/

        }/*FOR*/

      Free1D(cycles[i]);

    }/*FOR*/

  }/*FOR*/

  Free1D(counters);
  Free1D(path);
  Free1D(completed);
  Free1D(visited);
  Free1D(roots);

  for (int i = 0; i < ncycles; i++)
    Free1D(cycles[i]);
  Free1D(cycles);

  /* complete the path matrix by setting the diagonal. */
  for (int i = 0; i < nnodes; i++)
    final_rmat[CMC(i, i, nnodes)] = TRUE;

}/*COMPLETE_REACHABILITY_MATRIX*/

