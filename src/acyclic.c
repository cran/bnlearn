#include "common.h"

#define GOOD 1
#define BAD 0

static SEXP build_return_array(SEXP nodes, short int *status,
    int nrows, int check_status, SEXP return_nodes);

/*
 * Implementation  based on the proof of proposition 1.4.2
 * in "Digraphs Theory, Algorithms and Applications" by
 * Bang-Jensen and Gutin, page 13.
 */

SEXP is_pdag_acyclic(SEXP arcs, SEXP nodes, SEXP return_nodes, SEXP debug) {

int i = 0, j = 0, z = 0;
int nrows = LENGTH(nodes);
int check_status = nrows, check_status_old = nrows;
int *rowsums = NULL, *colsums = NULL, *crossprod = NULL, *a = NULL;
int *debuglevel = NULL;
short int *status = NULL;
SEXP amat;

  /* dereference the debug parameter. */
  debuglevel = LOGICAL(debug);

  /* build the adjacency matrix from the arc set.  */
  if (*debuglevel > 0)
    Rprintf("* building the adjacency matrix.\n");

  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* initialize the status, {row,col}sums and crossprod arrays. */
  status = allocstatus(nrows);
  rowsums = alloc1dcont(nrows);
  colsums = alloc1dcont(nrows);
  crossprod = alloc1dcont(nrows);

  if (*debuglevel > 0)
    Rprintf("* checking whether the partially directed graph is acyclic.\n");

  /* even in the worst case scenario at least two nodes are marked as
   * good at each iteration, so even ceil(nrows/2) iterations should be
   * enough. */
  for (z = 0; z < nrows; z++) {

start:

    if (*debuglevel > 0)
      Rprintf("* beginning iteration %d.\n", z + 1);

    for (i = 0; i < nrows; i++) {

      /* skip known-good nodes. */
      if (status[i] == GOOD) continue;

      /* reset and update row and column totals. */
      rowsums[i] = colsums[i] = crossprod[i] = 0;

      /* compute row and column totals for the i-th node. */
      for (j = 0; j < nrows; j++) {

        rowsums[i] += a[CMC(i, j, nrows)];
        colsums[i] += a[CMC(j, i, nrows)];
        crossprod[i] += a[CMC(i, j, nrows)] * a[CMC(j, i, nrows)];

      }/*FOR*/

there:

      if (*debuglevel > 0)
        Rprintf("  > checking node %s (%d child(ren), %d parent(s), %d neighbours).\n",
          NODE(i), rowsums[i], colsums[i], crossprod[i]);

      /* if either total is zero, the node is either a root node or a
       * leaf node, and is not part of any cycle. */
      if (((rowsums[i] == 0) || (colsums[i] == 0)) ||
          ((crossprod[i] == 1) && (rowsums[i] == 1) && (colsums[i] == 1))) {

        if (*debuglevel > 0)
          Rprintf("  @ node %s is cannot be part of a cycle.\n", NODE(i));

        /* update the adjacency matrix and the row/column totals. */
        for (j = 0; j < nrows; j++)
          a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;

        rowsums[i] = colsums[i] = crossprod[i] = 0;

        /* mark the node as good. */
        status[i] = GOOD;
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

          if (*debuglevel > 0)
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

      if (*debuglevel > 0)
        Rprintf("@ at least three nodes are needed to have a cycle.\n");

      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/

    /* if there are three or more bad nodes and there was no change in
     * the last iteration, the algorithm is stuck on a cycle. */
    if (check_status_old == check_status) {

      if (*debuglevel > 0)
        Rprintf("@ no change in the last iteration.\n");

      /* give up and call c_has_path() to kill some undirected arcs. */
      for (i = 0; i < nrows; i++)
        for (j = 0; j < i; j++)
          if (a[CMC(i, j, nrows)] * a[CMC(j, i, nrows)] == 1) {

            /* remove the arc from the adjacency matrix while testing it,
             * there's a path is always found (the arc itself). */
            a[CMC(i, j, nrows)] = a[CMC(j, i, nrows)] = 0;

            if(!c_has_path(i, j, INTEGER(amat), nrows, nodes, FALSE, TRUE, FALSE) &&
               !c_has_path(j, i, INTEGER(amat), nrows, nodes, FALSE, TRUE, FALSE)) {

              if (*debuglevel > 0)
                Rprintf("@ arc %s - %s is not part of any cycle, removing.\n", NODE(i), NODE(j));

              /* increase the iteration counter and start again. */
              z++;
              goto start;

            }/*THEN*/
            else {

              /* at least one cycle is really present; give up and return.  */
              UNPROTECT(1);
              return build_return_array(nodes, status, nrows, check_status, return_nodes);

            }/*ELSE*/

          }/*THEN*/

      /* give up if there are no undirected arcs, cycles composed
       * entirely by directed arcs are never false positives. */
      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/
    else {

      check_status_old = check_status;

    }/*ELSE*/

  }/*FOR*/

  UNPROTECT(1);
  return build_return_array(nodes, status, nrows, check_status, return_nodes);

}/*IS_PDAG_ACYCLIC*/

/* utility function to build the return value. */

static SEXP build_return_array(SEXP nodes, short int *status, int nrows,
    int check_status, SEXP return_nodes) {

int i = 0, j = 0;
SEXP res;

  if (check_status < 3) {

    if (isTRUE(return_nodes)) {

      PROTECT(res = allocVector(STRSXP, 0));

    }/*THEN*/
    else {

      PROTECT(res = allocVector(LGLSXP, 1));
      LOGICAL(res)[0] = TRUE;

    }/*ELSE*/

  }/*THEN*/
  else {

    if (isTRUE(return_nodes)) {

      PROTECT(res = allocVector(STRSXP, check_status));

      for (i = 0; i < nrows; i++)
        if (status[i] == BAD)
          SET_STRING_ELT(res, j++, STRING_ELT(nodes, i));

    }/*THEN*/
    else {

      PROTECT(res = allocVector(LGLSXP, 1));
      LOGICAL(res)[0] = FALSE;

    }/*ELSE*/

  }/*ELSE*/

  UNPROTECT(1);

  return res;

}/*BUILD_RETURN_ARRAY*/

