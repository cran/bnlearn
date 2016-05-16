#include "include/rcore.h"
#include "include/matrix.h"
#include "include/graph.h"

SEXP has_pdag_path(SEXP from, SEXP to, SEXP amat, SEXP nrows, SEXP nodes,
    SEXP underlying, SEXP exclude_direct, SEXP debug) {

int start = INT(from) - 1, stop = INT(to) - 1, n = INT(nrows);
int debuglevel = LOGICAL(debug)[0], notdirect = LOGICAL(exclude_direct)[0],
      ugraph = LOGICAL(underlying)[0], *a = INTEGER(amat);
int *path = NULL, *scratch = NULL, res = 0;

  /* allocate buffers for c_has_path(). */
  path = Calloc1D(n, sizeof(int));
  scratch = Calloc1D(n, sizeof(int));

  res = c_has_path(start, stop, a, n, nodes, ugraph, notdirect,
          path, scratch, debuglevel);

  Free1D(path);
  Free1D(scratch);

  return ScalarLogical(res);

}/*HAS_DAG_PATH*/

int c_has_path(int start, int stop, int *amat, int n, SEXP nodes,
    int ugraph, int notdirect, int *path, int *counter, int debuglevel) {

int i = 0, a1 = 0, a2 = 0, path_pos = 0, cur = start;

  /* remove any arc between start and stop if asked to. */
  if (notdirect) {

    a1 = amat[CMC(start, stop, n)];
    a2 = amat[CMC(stop, start, n)];
    amat[CMC(start, stop, n)] =  amat[CMC(stop, start, n)] = 0;

  }/*THEN*/

  /* initialize the position counters for the rows of the adjacency matrix. */
  memset(counter, '\0', n * sizeof(int));
  /* initialize the path array. */
  memset(path, '\0', n * sizeof(int));

  /* iterate until the other node is found. */
  while (cur != stop) {

    if (debuglevel > 0) {

      Rprintf("* currently at '%s'.\n", NODE(cur));
      Rprintf("  > current path is:\n");
      for (int i = 0; i < path_pos ; i ++)
        Rprintf("'%s' ", NODE(path[i]));
      Rprintf("'%s' \n", NODE(cur));

    }/*THEN*/

there:

    /* find the next child of the 'cur' node. */
    for (i = 0; (i < n) && (counter[cur] < n); i++) {

      if (!ugraph) {

        if (amat[CMC(cur, counter[cur], n)] != 0)
          break;

      }/*THEN*/
      else {

        /* if we are looking for a path in the underlying (undirected)
         * graph, check also the symmetric entry of the adjacency matrix. */
        if ((amat[CMC(cur, counter[cur], n)] != 0) ||
            (amat[CMC(counter[cur], cur, n)] != 0))
          break;

      }/*ELSE*/

      counter[cur]++;

    }/*FOR*/

    /* the column indexes range from 0 to n - 1;  counter value of n means
     * that all the children of that node has beeen visited. */
    if (counter[cur] == n) {

      /* if this node is the first one in the path, there search is finished
       * and the 'stop' node was not found; return FALSE. */
      if  (path_pos == 0) {

        /* remove any arc between start and stop if asked to. */
        if (notdirect) {

          amat[CMC(start, stop, n)] = a1;
          amat[CMC(stop, start, n)] = a2;

        }/*THEN*/

        return FALSE;

      }/*THEN*/

      if (debuglevel > 0)
        Rprintf("  > node '%s' has no more children, going back to '%s'.\n",
          NODE(cur), NODE(path[path_pos - 1]));

      /* this node has no more children, skip back to the previous one. */
      cur = path[--path_pos];
      path[path_pos + 1] = 0;

      goto there;

    }/*THEN*/
    else {

      /* inrement the counter to get to the next node to check, unless
       * that one was the last. */
      if (counter[cur] < n)
        counter[cur]++;

      /* do not visit an already visited node */
      for (i = path_pos - 1; i >= 0; i--) {

        if ((counter[cur] - 1) == path[i]) {

          if (debuglevel > 0)
            Rprintf("  @ node '%s' already visited, skipping.\n", NODE(path[i]));

          goto there;

        }/*THEN*/

      }/*FOR*/

      /* update the path. */
      path[path_pos++] = cur;
      /* the current node is now the children we have just found. */
      cur = counter[cur] - 1;

      if (debuglevel > 0)
        Rprintf("  > jumping to '%s'.\n", NODE(cur));

    }/*ELSE*/

  }/*WHILE*/


  /* remove any arc between start and stop if asked to. */
  if (notdirect) {

    amat[CMC(start, stop, n)] = a1;
    amat[CMC(stop, start, n)] = a2;

  }/*THEN*/

  /* node 'stop' has been found, return TRUE. */
  return TRUE;

}/*C_HAS_PATH*/

int c_directed_path(int start, int stop, int *amat, int n, SEXP nodes,
    int *path, int *counter, int debuglevel) {

int i = 0, path_pos = 0, cur = start;

  /* initialize the position counters for the rows of the adjacency matrix. */
  memset(counter, '\0', n * sizeof(int));
  /* initialize the path array. */
  memset(path, '\0',  n * sizeof(int));

  /* iterate until the other node is found. */
  while (cur != stop) {

    if (debuglevel > 0) {

      Rprintf("* currently at '%s'.\n", NODE(cur));
      Rprintf("  > current path is:\n");
      for (int i = 0; i < path_pos ; i ++)
        Rprintf("'%s' ", NODE(path[i]));
      Rprintf("'%s' \n", NODE(cur));

    }/*THEN*/

there:

    /* find the next child of the 'cur' node. */
    for (i = 0; (i < n) && (counter[cur] < n); i++) {

      if ((amat[CMC(cur, counter[cur], n)] != 0) &&
          (amat[CMC(counter[cur], cur, n)] == 0))
        break;

      counter[cur]++;

    }/*FOR*/

    /* the column indexes range from 0 to n - 1;  counter value of n means
     * that all the children of that node has beeen visited. */
    if (counter[cur] == n) {

      /* if this node is the first one in the path, there search is finished
       * and the 'stop' node was not found; return FALSE. */
      if  (path_pos == 0)
        return FALSE;

      if (debuglevel > 0)
        Rprintf("  > node '%s' has no more children, going back to '%s'.\n",
          NODE(cur), NODE(path[path_pos - 1]));

      /* this node has no more children, skip back to the previous one. */
      cur = path[--path_pos];
      path[path_pos + 1] = 0;

      goto there;

    }/*THEN*/
    else {

      /* inrement the counter to get to the next node to check, unless
       * that one was the last. */
      if (counter[cur] < n)
        counter[cur]++;

      /* do not visit an already visited node */
      for (i = path_pos - 1; i >= 0; i--) {

        if ((counter[cur] - 1) == path[i]) {

          if (debuglevel > 0)
            Rprintf("  @ node '%s' already visited, skipping.\n", NODE(path[i]));

          goto there;

        }/*THEN*/

      }/*FOR*/

      /* update the path. */
      path[path_pos++] = cur;
      /* the current node is now the children we have just found. */
      cur = counter[cur] - 1;

      if (debuglevel > 0)
        Rprintf("  > jumping to '%s'.\n", NODE(cur));

    }/*ELSE*/

  }/*WHILE*/

  /* node 'stop' has been found, return TRUE. */
  return TRUE;

}/*C_DIRECTED_PATH*/

int c_uptri3_path(short int *uptri, int *depth, int from, int to, int nnodes,
    SEXP nodes, int debuglevel) {

int i = 0, j = 0, d = 0;

  /* reset the depth counter for all the nodes. */
  memset(depth, '\0', nnodes * sizeof(int));
  depth[from] = 1;

  /* for each depth level... */
  for (d = 1; d <= nnodes; d++) {

    if (debuglevel)
      Rprintf("* considering depth %d.\n", d);

    /* ... for each node... */
    for (i = 0; i < nnodes; i++) {

      /* ... if it is at the current depth level... */
      if (depth[i] != d)
        continue;

      if (debuglevel)
        Rprintf("  > found node %s.\n", NODE(i));

      for (j = 0; j < nnodes; j++) {

        /* ... not on the diagonal (which is not allocated by UPTRI) ... */
        if (i == j)
          continue;
        /* ... and it is adjacent to another node... */
        if (!(uptri[UPTRI3(i + 1, j + 1, nnodes)] == 1))
          continue;

        /* ... that hasn't already been visited... */
        if (depth[j] != 0) {

          if (debuglevel)
            Rprintf("  @ node '%s' already visited, skipping.\n", NODE(j));

          continue;

        }/*THEN*/

        if (j == to) {

          /* ... and it's the destination, exit. */
          if (debuglevel)
            Rprintf("  @ arrived at node %s, exiting.\n", NODE(to));

          return TRUE;

        }/*THEN*/
        else {

          /* ...and it's not the destination, add it to the next depth level. */
          depth[j] = d + 1;

        }/*ELSE*/

        if (debuglevel)
          Rprintf("  > added node %s at depth %d\n", NODE(j), d + 1);

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  return FALSE;

}/*C_UPTRI3_PATH*/

