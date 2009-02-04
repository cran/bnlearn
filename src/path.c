#include "common.h"

#define AMAT(i,j) INTEGER(amat)[i + j * n]
#define NODE(i) CHAR(STRING_ELT(nodes, i))

SEXP has_dag_path(SEXP from, SEXP to, SEXP amat, SEXP nrows, SEXP nodes, 
    SEXP underlying, SEXP debug) {

  int start = INTEGER(from)[0] - 1;
  int stop = INTEGER(to)[0] - 1;
  int n = INTEGER(nrows)[0];
  int *counter;
  int *path;
  int path_pos = 0;
  int cur = start;
  int i = 0;

  SEXP res;

  /* initialize the position counters for the rows of the adjacency matrix. */
  counter = (int *) R_alloc(n, sizeof(int));
  memset(counter, '\0', sizeof(int) * n);
  /* initialize the path array. */
  path = (int *) R_alloc(n, sizeof(int));
  memset(path, '\0', sizeof(int) * n);

  /* allocate the result. */
  PROTECT(res = allocVector(LGLSXP, 1));
  LOGICAL(res)[0] = TRUE;

  /* iterate until the other node is found. */
  while (cur != stop) {

    if (isTRUE(debug)) {

      Rprintf("* currently at '%s'.\n", NODE(cur));
      Rprintf("  > current path is:\n");
      for (int i = 0; i < path_pos ; i ++)
        Rprintf("'%s' ", NODE(path[i]));
      Rprintf("'%s' \n", NODE(cur));

    }/*THEN*/

there:

    /* find the next child of the 'cur' node. */
    for (i = 0; (i < n) && (counter[cur] < n); i++) {

      if (!isTRUE(underlying)) {

        if (AMAT(cur, counter[cur]) != 0)
          break;

      }/*THEN*/
      else {

        /* if we are looking for a path in the underlying (undirected)
         * graph, check also the symmetric entry of the adjacency matrix. */
        if ((AMAT(cur, counter[cur]) != 0) || (AMAT(counter[cur], cur) != 0))
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

        LOGICAL(res)[0] = FALSE;
        UNPROTECT(1);
        return res;

      }/*THEN*/

      if (isTRUE(debug)) {

        Rprintf("  > node '%s' has no more children, going back to '%s'.\n",
          NODE(cur), NODE(path[path_pos - 1]));

      }/*THEN*/

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

          if (isTRUE(debug)) {

            Rprintf("  @ node '%s' already visited, skipping.\n", NODE(path[i]));

          }/*THEN*/

          goto there;

        }/*THEN*/

      }/*FOR*/

      /* update the path. */
      path[path_pos++] = cur;
      /* the current node is now the children we have just found. */
      cur = counter[cur] - 1;

      if (isTRUE(debug)) {

        Rprintf("  > jumping to '%s'.\n", NODE(cur));

      }/*THEN*/

    }/*ELSE*/

  }/*WHILE*/

  /* node 'stop' has been found, return TRUE. */
  UNPROTECT(1);
  return res;

}/*HAS_DAG_PATH*/

SEXP how_many_cycles(SEXP from, SEXP to, SEXP amat, SEXP nrows, SEXP nodes, SEXP debug) {

  int start = INTEGER(from)[0] - 1;
  int stop = INTEGER(to)[0] - 1;
  int n = INTEGER(nrows)[0];
  int *counter;
  int *path;
  int path_pos = 0;
  int cur = start;
  int cycle_counter = 0;
  int i = 0;

  SEXP res;

  /* initialize the position counters for the rows of the adjacency matrix. */
  counter = (int *) R_alloc(n, sizeof(int));
  memset(counter, '\0', sizeof(int) * n);
  /* initialize the path array. */
  path = (int *) R_alloc(n, sizeof(int));
  memset(path, '\0', sizeof(int) * n);

  /* allocate the result. */
  PROTECT(res = allocVector(INTSXP, 1));
  INTEGER(res)[0] = 0;

  while (1) {

    if (isTRUE(debug))
      Rprintf("* currently at '%s'.\n", NODE(cur));

there:

    /* find the next child of the 'cur' node. */
    for (i = 0; (i < n) && (counter[cur] < n); i++) {

      if (AMAT(cur, counter[cur]) != 0)
        break;

      counter[cur]++;

    }/*FOR*/

    if (cur == stop) {

      cycle_counter++;

      if (isTRUE(debug)) {

        Rprintf("  @ found node '%s' ! cycle counter is now %d.\n",
          NODE(stop), cycle_counter);
        Rprintf("  @ the cycle is:\n'%s' ", NODE(stop));
        for (int i = 0; i < path_pos ; i ++)
          Rprintf("'%s' ", NODE(path[i]));
        Rprintf("'%s' \n", NODE(cur));

      }/*THEN*/

      /* do not go on on this path. */
      if (counter[cur] < n)
        counter[cur]++;
      cur = path[--path_pos];
      path[path_pos + 1] = 0;

    }/*THEN*/
    else if (counter[cur] == n) {

      /* the column indexes range from 0 to n - 1;  counter value of n means
       * that all the children of that node has beeen visited.  */

      /* if this node is the first one in the path, there search is finished
       * and the 'stop' node was not found; return FALSE.*/
      if  (path_pos == 0) {

        INTEGER(res)[0] = cycle_counter;
        UNPROTECT(1);
        return res;

      }/*THEN*/

      if (isTRUE(debug)) {

        Rprintf("  > node '%s' has no more children, going back to '%s'.\n",
          NODE(cur), NODE(path[path_pos - 1]));

      }/*THEN*/

      /* this node has no more children, skip back to the previous one. */
      cur = path[--path_pos];
      path[path_pos + 1] = 0;

    }/*THEN*/
    else {

      /* inrement the counter to get to the next node to check, unless
       * that one was the last. */
      if (counter[cur] < n)
        counter[cur]++;

      /* do not visit an already visited node */
      for (i = path_pos - 1; i >= 0; i--) {

        if ((counter[cur] - 1) == path[i]) {

          if (isTRUE(debug)) {

            Rprintf("  @ node '%s' already visited, skipping.\n", NODE(path[i]));

          }/*THEN*/

          goto there;

        }/*THEN*/

      }/*FOR*/

      /* update the path. */
      path[path_pos++] = cur;
      /* the current node is now the children we have just found. */
      cur = counter[cur] - 1;

      if (isTRUE(debug)) {

        Rprintf("  > jumping to '%s'.\n", NODE(cur));

      }/*THEN*/

    }/*ELSE*/

  }/*WHILE*/

}/*HOW_MANY_CYCLES*/

