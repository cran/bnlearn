#include "include/rcore.h"
#include "include/matrix.h"
#include "include/graph.h"

/* utility function to build the return value. */
static SEXP build_return_array(SEXP nodes, short int *skip, int nrows,
    int check_status, SEXP return_nodes) {

int i = 0, j = 0;
SEXP res;

  if (check_status < 3) {

    if (isTRUE(return_nodes))
      return allocVector(STRSXP, 0);
    else
      return ScalarLogical(TRUE);

  }/*THEN*/
  else {

    if (isTRUE(return_nodes)) {

      PROTECT(res = allocVector(STRSXP, check_status));

      for (i = 0; i < nrows; i++)
        if (skip[i] == FALSE)
          SET_STRING_ELT(res, j++, STRING_ELT(nodes, i));

      UNPROTECT(1);

      return res;

    }/*THEN*/
    else {

      return ScalarLogical(FALSE);

    }/*ELSE*/

  }/*ELSE*/

}/*BUILD_RETURN_ARRAY*/

/* implementation based on the proof of proposition 1.4.2 in "Digraphs Theory, 
 * Algorithms and Applications" by Bang-Jensen and Gutin, page 13. */
SEXP is_pdag_acyclic(SEXP arcs, SEXP nodes, SEXP return_nodes,
    SEXP directed, SEXP debug) {

int i = 0, j = 0, z = 0;
int nrows = length(nodes);
int check_status = nrows, check_status_old = nrows;
int *rowsums = NULL, *colsums = NULL, *crossprod = NULL, *a = NULL;
int *path = NULL, *scratch = NULL, debuglevel = isTRUE(debug);
short int *skip = NULL;
SEXP amat, res;

  /* build the adjacency matrix from the arc set.  */
  if (debuglevel > 0)
    Rprintf("* building the adjacency matrix.\n");

  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* should we consider only directed arcs? */
  if (isTRUE(directed)) {

    /* removing undirected arcs, so that only cycles made only of directed
     * arcs will make the function return TRUE. */

    for (i = 0; i < nrows; i++)
      for (j = 0; j < nrows; j++)
        if ((a[CMC(i, j, nrows)] == 1) && (a[CMC(j, i, nrows)] == 1))
          a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;

  }/*THEN*/

  /* initialize the status, {row,col}sums and crossprod arrays. */
  skip = Calloc1D(nrows, sizeof(short int));
  rowsums = Calloc1D(nrows, sizeof(int));
  colsums = Calloc1D(nrows, sizeof(int));
  crossprod = Calloc1D(nrows, sizeof(int));

  /* allocate buffers for c_has_path(). */
  path = Calloc1D(nrows, sizeof(int));
  scratch = Calloc1D(nrows, sizeof(int));

  if (debuglevel > 0)
    Rprintf("* checking whether the partially directed graph is acyclic.\n");

  /* even in the worst case scenario at least two nodes are marked as
   * good at each iteration, so even ceil(nrows/2) iterations should be
   * enough. */
  for (z = 0; z < nrows; z++) {

start:

    if (debuglevel > 0)
      Rprintf("* beginning iteration %d.\n", z + 1);

    for (i = 0; i < nrows; i++) {

      /* skip known-good nodes. */
      if (skip[i] == TRUE) continue;

      /* reset and update row and column totals. */
      rowsums[i] = colsums[i] = crossprod[i] = 0;

      /* compute row and column totals for the i-th node. */
      for (j = 0; j < nrows; j++) {

        rowsums[i] += a[CMC(i, j, nrows)];
        colsums[i] += a[CMC(j, i, nrows)];
        crossprod[i] += a[CMC(i, j, nrows)] * a[CMC(j, i, nrows)];

      }/*FOR*/

there:

      if (debuglevel > 0)
        Rprintf("  > checking node %s (%d child(ren), %d parent(s), %d neighbours).\n",
          NODE(i), rowsums[i], colsums[i], crossprod[i]);

      /* if either total is zero, the node is either a root node or a
       * leaf node, and is not part of any cycle. */
      if (((rowsums[i] == 0) || (colsums[i] == 0)) ||
          ((crossprod[i] == 1) && (rowsums[i] == 1) && (colsums[i] == 1))) {

        if (debuglevel > 0)
          Rprintf("  @ node %s is cannot be part of a cycle.\n", NODE(i));

        /* update the adjacency matrix and the row/column totals. */
        for (j = 0; j < nrows; j++)
          a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;

        rowsums[i] = colsums[i] = crossprod[i] = 0;

        /* mark the node as good. */
        skip[i] = TRUE;
        check_status--;

      }/*THEN*/
      else if (crossprod[i] == 1) {

        /* find the other of the undirected arc. */
        for (j = 0; j < i; j++)
          if (a[CMC(i, j, nrows)] * a[CMC(j, i, nrows)] == 1)
            break;

        /* safety check, just in case. */
        if (i == j) continue;

        if (((colsums[i] == 1) && (colsums[j] == 1)) ||
            ((rowsums[i] == 1) && (rowsums[j] == 1))) {

          if (debuglevel > 0)
            Rprintf("  @ arc %s - %s is cannot be part of a cycle.\n", NODE(i), NODE(j));

          /* update the adjacency matrix and the row/column totals. */
          a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;
          crossprod[i] = 0;
          rowsums[i]--;
          colsums[i]--;
          rowsums[j]--;
          colsums[j]--;

          /* jump back to the first check; if either the row or column total
           * was equal to 1 only because of the undirected arc, the node can
           * now be marked as good. */
          if ((rowsums[i] == 0) || (colsums[i] == 0))
            goto there;

        }/*THEN*/

      }/*THEN*/

    }/*FOR*/

    /* at least three nodes are needed to have a cycle. */
    if (check_status < 3) {

      if (debuglevel > 0)
        Rprintf("@ at least three nodes are needed to have a cycle.\n");

      goto end;

    }/*THEN*/

    /* if there are three or more bad nodes and there was no change in
     * the last iteration, the algorithm is stuck on a cycle. */
    if (check_status_old == check_status) {

      if (debuglevel > 0)
        Rprintf("@ no change in the last iteration.\n");

      /* give up and call c_has_path() to kill some undirected arcs. */
      for (i = 0; i < nrows; i++)
        for (j = 0; j < i; j++)
          if (a[CMC(i, j, nrows)] * a[CMC(j, i, nrows)] == 1) {

            /* remove the arc from the adjacency matrix while testing it,
             * there's a path is always found (the arc itself). */
            a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;

            if(!c_has_path(i, j, INTEGER(amat), nrows, nodes, FALSE, TRUE, path, scratch, FALSE) &&
               !c_has_path(j, i, INTEGER(amat), nrows, nodes, FALSE, TRUE, path, scratch, FALSE)) {

              if (debuglevel > 0)
                Rprintf("@ arc %s - %s is not part of any cycle, removing.\n", NODE(i), NODE(j));

              /* increase the iteration counter and start again. */
              z++;
              goto start;

            }/*THEN*/
            else {

              /* at least one cycle is really present; give up and return.  */
              goto end;

            }/*ELSE*/

          }/*THEN*/

      /* give up if there are no undirected arcs, cycles composed entirely by
       * directed arcs are never false positives. */
      goto end;

    }/*THEN*/
    else {

      check_status_old = check_status;

    }/*ELSE*/

  }/*FOR*/

end:

  res = build_return_array(nodes, skip, nrows, check_status, return_nodes);

  Free1D(skip);
  Free1D(rowsums);
  Free1D(colsums);
  Free1D(crossprod);
  Free1D(path);
  Free1D(scratch);

  UNPROTECT(1);

  return res;

}/*IS_PDAG_ACYCLIC*/

