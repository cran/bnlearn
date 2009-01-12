#include "common.h"

#define NODE(i) CHAR(STRING_ELT(nodes, i))
#define isTRUE(logical) LOGICAL(logical)[0] == TRUE
#define AMAT(i,j) INTEGER(amat)[i + j * nrows]

#define BLANKET		 1
#define NEIGHBOUR	 2
#define PARENT		 3
#define CHILD		 4

SEXP cache_node_structure(int cur, SEXP nodes, SEXP amat, int nrows, SEXP debug) {

  int i = 0, j = 0;
  SEXP structure, names;
  SEXP mb, nbr, children, parents;
  int num_parents = 0, num_children = 0, num_neighbours = 0, num_blanket = 0;
  unsigned short int *status;

  /* allocate and initialize the status vector. */
  status = (unsigned short int *) R_alloc(nrows, sizeof(short int));
  memset(status, '\0', sizeof(short int) * nrows);

  if (isTRUE(debug))
    Rprintf("* node %s.\n", NODE(cur));

  for (i = 0; i < nrows; i++) {

    if (AMAT(cur, i) == 1) {

      if (AMAT(i, cur) == 0) {

        /* if a[i,j] = 1 and a[j,i] = 0, then i -> j. */
        if (isTRUE(debug))
          Rprintf("  > found child %s.\n", NODE(i));

        num_children++;
        num_neighbours++;
        status[i] = CHILD;

        /* check whether this child has any other parent. */
        for (j = 0; j < nrows; j++) {

          if ((AMAT(j,i) == 1) && (AMAT(i,j) == 0) && (j != cur)) {

            if (isTRUE(debug))
              Rprintf("  > found node %s in markov blanket.\n", NODE(j));

            num_blanket++;
            status[j] = BLANKET;

          }/*THEN*/

        }/*FOR*/       

      }/*THEN*/
      else {

        /* if a[i,j] = 1 and a[j,i] = 1, then i -- j. */
        if (isTRUE(debug))
          Rprintf("  > found neighbour %s.\n", NODE(i));

        num_neighbours++;
        status[i] = NEIGHBOUR;

      }/*ELSE*/

    }/*THEN*/
    else {

      if (AMAT(i, cur) == 1) {

        /* if a[i,j] = 0 and a[j,i] = 1, then i <- j. */
        if (isTRUE(debug))
          Rprintf("  > found parent %s.\n", NODE(i));

        num_parents++;
        num_neighbours++;
        status[i] = PARENT;

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  if (isTRUE(debug))
    Rprintf("  > node %s has %d parent(s), %d child(ren), %d neighbour(s) and %d nodes in the markov blanket.\n", 
      NODE(cur), num_parents, num_children, num_neighbours, num_blanket + num_neighbours);

  /* allocate and initialize the names of the elements. */
  PROTECT(names = allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, mkChar("mb"));
  SET_STRING_ELT(names, 1, mkChar("nbr"));
  SET_STRING_ELT(names, 2, mkChar("parents"));
  SET_STRING_ELT(names, 3, mkChar("children"));

  /* allocate the list and set its attributes. */
  PROTECT(structure = allocVector(VECSXP, 4));
  setAttrib(structure, R_NamesSymbol, names);

  /* allocate and fill the "children" element of the list. */
  PROTECT(children = allocVector(STRSXP, num_children));
  for (i = 0, j = 0; (i < nrows) && (j < num_children); i++) {

    if (status[i] == CHILD)
      SET_STRING_ELT(children, j++, STRING_ELT(nodes, i));

  }/*FOR*/

  /* allocate and fill the "parents" element of the list. */
  PROTECT(parents = allocVector(STRSXP, num_parents));
  for (i = 0, j = 0; (i < nrows) && (j < num_parents); i++) {

    if (status[i] == PARENT)
      SET_STRING_ELT(parents, j++, STRING_ELT(nodes, i));

  }/*FOR*/

  /* allocate and fill the "nbr" element of the list. */
  PROTECT(nbr = allocVector(STRSXP, num_neighbours));
  for (i = 0, j = 0; (i < nrows) && (j < num_neighbours); i++) {

    if (status[i] >= NEIGHBOUR)
      SET_STRING_ELT(nbr, j++, STRING_ELT(nodes, i));

  }/*FOR*/

  /* allocate and fill the "mb" element of the list. */
  PROTECT(mb = allocVector(STRSXP, num_blanket + num_neighbours));
  for (i = 0, j = 0; (i < nrows) && (j < num_blanket + num_neighbours); i++) {

    if (status[i] >= BLANKET)
      SET_STRING_ELT(mb, j++, STRING_ELT(nodes, i));

  }/*FOR*/

  /* attach the string vectors to the list. */
  SET_VECTOR_ELT(structure, 0, mb); 
  SET_VECTOR_ELT(structure, 1, nbr); 
  SET_VECTOR_ELT(structure, 2, parents); 
  SET_VECTOR_ELT(structure, 3, children); 

  UNPROTECT(6);

  return structure;

}/*CACHE_NODE_STRUCTURE*/

SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug) {

  int i = 0;
  int length_nodes = LENGTH(nodes);

  SEXP bn, temp;

  /* allocate the list and set its attributes.*/
  PROTECT(bn = allocVector(VECSXP, length_nodes));
  setAttrib(bn, R_NamesSymbol, nodes);

  if (isTRUE(debug))
    Rprintf("* (re)building cached information about network structure.\n");

  /* populate the list with nodes' data. */
  for (i = 0; i < length_nodes; i++) {

    temp = cache_node_structure(i, nodes, amat, length_nodes, debug);

    /* save the returned list. */
    SET_VECTOR_ELT(bn, i, temp);

  }/*FOR*/

  UNPROTECT(1);

  return bn;

}/*CACHE_STRUCTURE*/

