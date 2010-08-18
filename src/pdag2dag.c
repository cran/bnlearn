#include "common.h"

/* return the complete orientation of a graph (the nodes argument gives
  * the node ordering). */
SEXP pdag2dag(SEXP arcs, SEXP nodes) {

int i = 0, j = 0, n = LENGTH(nodes);
int *a = NULL;
SEXP amat, res;

  /* build the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* scan the adjacency matrix. */
  for (i = 0; i < n; i++) {

    for (j = i + 1; j < n; j++) {

      /* if an arc is undirected, kill the orientation that violates the
       * specified node ordering (the one which is located in the lower
       * half of the matrix). */
      if ((a[CMC(i, j, n)] == 1) && (a[CMC(j, i, n)] == 1))
        a[CMC(j, i, n)] = 0;

    }/*FOR*/

  }/*FOR*/

  /* build the return value. */
  PROTECT(res = amat2arcs(amat, nodes));

  UNPROTECT(2);

  return res;

}/*PDAG2DAG*/

/* return the skeleton of a DAG/PDAG. */
SEXP dag2ug(SEXP bn, SEXP moral, SEXP debug) {

int i = 0, j = 0, k = 0, nnodes = 0, narcs = 0, row = 0;
int *debuglevel = LOGICAL(debug), *moralize = LOGICAL(moral);\
int *nparents = NULL, *nnbr = NULL;
SEXP node_data, current, nodes, result, temp;

  /* get the nodes' data. */
  node_data = getListElement(bn, "nodes");
  nnodes = LENGTH(node_data);
  nodes = getAttrib(node_data, R_NamesSymbol);

  /* allocate and initialize parents' and neighbours' counters. */
  nnbr = alloc1dcont(nnodes);
  if (*moralize > 0)
    nparents = alloc1dcont(nnodes);

  /* first pass: count neighbours, parents and resulting arcs. */
  for (i = 0; i < nnodes; i++) {

    /* get the number of neighbours.  */
    current = VECTOR_ELT(node_data, i);
    nnbr[i] = LENGTH(getListElement(current, "nbr"));

    /* update the number of arcs to be returned. */
    if (*moralize > 0) {

      /* get also the number of parents, needed to account for the arcs added
       * for their moraliztion. */
      nparents[i] = LENGTH(getListElement(current, "parents"));
      narcs += nnbr[i] + nparents[i] * (nparents[i] - 1);

    }/*THEN*/
    else {

      narcs += nnbr[i];

    }/*ELSE*/

    if (*debuglevel > 0)  {

      if (*moralize > 0) {

        Rprintf("* scanning node %s, found %d neighbours and %d parents.\n",
          NODE(i), nnbr[i], nparents[i]);
        Rprintf("  > adding %d arcs, for a total of %d.\n",
          nnbr[i] + nparents[i] * (nparents[i] - 1), narcs);

      }/*THEN*/
      else {

        Rprintf("* scanning node %s, found %d neighbours.\n", NODE(i), nnbr[i]);
        Rprintf("  > adding %d arcs, for a total of %d.\n", nnbr[i], narcs);

      }/*ELSE*/

    }/*THEN*/

  }/*FOR*/

  /* allocate the return value. */
  PROTECT(result = allocMatrix(STRSXP, narcs, 2));
  /* allocate and set the column names. */
  finalize_arcs(result);

  /* second pass: fill the return value. */
  for (i = 0; i < nnodes; i++) {

    /* get to the current node. */
    current = VECTOR_ELT(node_data, i);
    /* get the neighbours. */
    temp = getListElement(current, "nbr");

    for (j = 0; j < nnbr[i]; j++) {

      SET_STRING_ELT(result, CMC(row, 0, narcs), STRING_ELT(nodes, i));
      SET_STRING_ELT(result, CMC(row, 1, narcs), STRING_ELT(temp, j));
      row++;

    }/*FOR*/

    /* if we are not creating a moral graph we are done with this node. */
    if (*moralize == 0)
      continue;

    /* get the parents. */
    temp = getListElement(current, "parents");

    for (j = 0; j < nparents[i]; j++) {

      for (k = j+1; k < nparents[i]; k++) {

        SET_STRING_ELT(result, CMC(row, 0, narcs), STRING_ELT(temp, k));
        SET_STRING_ELT(result, CMC(row, 1, narcs), STRING_ELT(temp, j));
        row++;
        SET_STRING_ELT(result, CMC(row, 0, narcs), STRING_ELT(temp, j));
        SET_STRING_ELT(result, CMC(row, 1, narcs), STRING_ELT(temp, k));
        row++;

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);

  /* be really sure not to return duplicate arcs in moral graphs when shielded
   * parents are present (the "shielding" are is counted twice). */
  if (*moralize > 0)
    return c_unique_arcs(result, nodes, FALSE);
  else
    return result;

}/*DAG2UG*/

