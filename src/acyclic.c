#include "common.h"

#define NODE(i) CHAR(STRING_ELT(nodes, i))
#define AMAT(i,j) INTEGER(amat)[i + j * nrows]
#define GOOD 1
#define BAD 0

static SEXP build_return_array(SEXP nodes, short int *status,
    int nrows, int check_status, SEXP return_nodes);

/*
 * Implementation  based on the proof of proposition 1.4.2
 * in "Digraphs Theory, Algorithms and Applications" by
 * Bang-Jensen and Gutin, page 13.
 */

SEXP is_dag_acyclic(SEXP arcs, SEXP nodes, SEXP return_nodes, SEXP debug) {

  SEXP amat;
  int i = 0, j = 0, z = 0;
  int nrows = LENGTH(nodes);
  int check_status = nrows;
  int check_status_old = nrows;
  int rowsums, colsums;
  short int *status;

  /* build the adjacency matrix from the arc set.  */
  if (isTRUE(debug))
    Rprintf("* building the adjacency matrix.\n");

  PROTECT(amat = arcs2amat(arcs, nodes));

  /* allocate and initialize the status array. */
  status = allocstatus(nrows);

  if (isTRUE(debug))
    Rprintf("* checking whether the directed graph is acyclic.\n");

  /* even in the worst case scenario at least two nodes are marked as
   * good at each iteration, so even ceil(nrows/2) iterations should be
   * enough. */
  for (z = 0; z < nrows; z++) {

    if (isTRUE(debug))
      Rprintf("* beginning iteration %d.\n", z + 1);

    for (i = 0; i < nrows; i++) {

      /* skip known-good nodes. */
      if (status[i] == GOOD) continue;

      /* reset and update row and column totals. */
      rowsums = colsums = 0;

      /* compute row and column totals for the i-th node. */
      for (j = 0; j < nrows; j++) {

        rowsums += AMAT(i, j);
        colsums += AMAT(j, i);

      }/*FOR*/

      if (isTRUE(debug))
        Rprintf("  > checking node %s (%d child(ren), %d parent(s)).\n",
          NODE(i), rowsums, colsums);

      /* if either total is zero, the node is either a root node or a
       * leaf node, and is not part of any cycle. */
      if ((rowsums == 0) || (colsums == 0)) {

        if (isTRUE(debug))
          Rprintf("  @ node %s is cannot be part of a cycle.\n", NODE(i));

        for (j = 0; j < nrows; j++)
          AMAT(i, j) = AMAT(j, i) = 0;

        /* mark the node as good. */
        status[i] = GOOD;
        check_status--;

      }/*THEN*/

    }/*FOR*/

    /* at least three nodes are needed to have a cycle. */
    if (check_status < 3) {

      if (isTRUE(debug))
        Rprintf("@ at least three nodes are needed to have a cycle.\n");

      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/

    /* if there are three or more bad nodes and there was no change in
     * the last iteration, the algorithm is stuck on a cycle. */
    if (check_status_old == check_status) {

      if (isTRUE(debug))
        Rprintf("@ no change in the last iteration.\n");

      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/
    else
      check_status_old = check_status;

  }/*FOR*/

  UNPROTECT(1);
  return build_return_array(nodes, status, nrows, check_status, return_nodes);

}/*IS_DAG_ACYCLIC*/

/*
 * An extension of the previous algorithm to the mixed graph case.
 */

SEXP is_pdag_acyclic(SEXP arcs, SEXP nodes, SEXP return_nodes, SEXP debug) {

  SEXP amat;
  int i = 0, j = 0, z = 0;
  int nrows = LENGTH(nodes);
  int check_status = nrows;
  int check_status_old = nrows;
  int *rowsums, *colsums, *crossprod;
  short int *status;

  /* build the adjacency matrix from the arc set.  */
  if (isTRUE(debug))
    Rprintf("* building the adjacency matrix.\n");

  PROTECT(amat = arcs2amat(arcs, nodes));

  /* initialize the status, {row,col}sums and crossprod arrays. */
  status = allocstatus(nrows);
  rowsums = alloc1dcont(nrows);
  colsums = alloc1dcont(nrows);
  crossprod = alloc1dcont(nrows);

  if (isTRUE(debug))
    Rprintf("* checking whether the partially directed graph is acyclic.\n");

  /* even in the worst case scenario at least two nodes are marked as
   * good at each iteration, so even ceil(nrows/2) iterations should be
   * enough. */
  for (z = 0; z < nrows; z++) {

    if (isTRUE(debug))
      Rprintf("* beginning iteration %d.\n", z + 1);

    for (i = 0; i < nrows; i++) {

      /* skip known-good nodes. */
      if (status[i] == GOOD) continue;

      /* reset and update row and column totals. */
      rowsums[i] = colsums[i] = crossprod[i] = 0;

      /* compute row and column totals for the i-th node. */
      for (j = 0; j < nrows; j++) {

        rowsums[i] += AMAT(i, j);
        colsums[i] += AMAT(j, i);
        crossprod[i] += AMAT(i, j) * AMAT(j, i);

      }/*FOR*/

there:

      if (isTRUE(debug))
        Rprintf("  > checking node %s (%d child(ren), %d parent(s), %d neighbours).\n",
          NODE(i), rowsums[i], colsums[i], crossprod[i]);

      /* if either total is zero, the node is either a root node or a
       * leaf node, and is not part of any cycle. */
      if (((rowsums[i] == 0) || (colsums[i] == 0)) ||
          ((crossprod[i] == 1) && (rowsums[i] == 1) && (colsums[i] == 1))) {

        if (isTRUE(debug))
          Rprintf("  @ node %s is cannot be part of a cycle.\n", NODE(i));

        /* update the adjacency matrix and the row/column totals. */
        for (j = 0; j < nrows; j++)
          AMAT(i, j) = AMAT(j, i) = 0;

        rowsums[i] = colsums[i] = crossprod[i] = 0;

        /* mark the node as good. */
        status[i] = GOOD;
        check_status--;

      }/*THEN*/
      else if (crossprod[i] == 1) {

        /* find the other of the undirected arc. */
        for (j = 0; j < i; j++)
          if (AMAT(i, j) * AMAT(j, i) == 1)
            break;

        /* safety check, just in case. */
        if (i == j) continue;

        if (((colsums[i] == 1) && (colsums[j] == 1)) ||
            ((rowsums[i] == 1) && (rowsums[j] == 1))) {

          if (isTRUE(debug))
            Rprintf("  @ arc %s - %s is cannot be part of a cycle.\n", NODE(i), NODE(j));

          /* update the adjacency matrix and the row/column totals. */
          AMAT(i, j) = AMAT(j, i) = 0;
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

      if (isTRUE(debug))
        Rprintf("@ at least three nodes are needed to have a cycle.\n");

      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/

    /* if there are three or more bad nodes and there was no change in
     * the last iteration, the algorithm is stuck on a cycle. */
    if (check_status_old == check_status) {

      if (isTRUE(debug))
        Rprintf("@ no change in the last iteration.\n");

      UNPROTECT(1);
      return build_return_array(nodes, status, nrows, check_status, return_nodes);

    }/*THEN*/
    else
      check_status_old = check_status;

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

