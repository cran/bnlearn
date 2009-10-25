#include "common.h"

#define BLANKET		 1
#define NEIGHBOUR	 2
#define PARENT		 3
#define CHILD		 4

static SEXP cache_node_structure(int cur, SEXP nodes, int *amat, int nrows,
    int *status, int debuglevel);

/* compute the cached values fro all nodes. */
SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug) {

int i = 0, debuglevel = LOGICAL(debug)[0], length_nodes = LENGTH(nodes);
int *status = NULL, *a = INTEGER(amat);

  SEXP bn, temp;

  /* allocate the list and set its attributes.*/
  PROTECT(bn = allocVector(VECSXP, length_nodes));
  setAttrib(bn, R_NamesSymbol, nodes);

  /* allocate and intialize the status vector. */
  status = alloc1dcont(length_nodes);

  if (isTRUE(debug))
    Rprintf("* (re)building cached information about network structure.\n");

  /* populate the list with nodes' data. */
  for (i = 0; i < length_nodes; i++) {

    /* (re)initialize the status vector. */
    memset(status, '\0', sizeof(int) * length_nodes);

    temp = cache_node_structure(i, nodes, a, length_nodes, status, debuglevel);

    /* save the returned list. */
    SET_VECTOR_ELT(bn, i, temp);

  }/*FOR*/

  UNPROTECT(1);

  return bn;

}/*CACHE_STRUCTURE*/

/* compute the cached values for a single node (R-friendly). */
SEXP cache_partial_structure(SEXP nodes, SEXP target, SEXP amat, SEXP debug) {

int i = 0, debuglevel = LOGICAL(debug)[0], length_nodes = LENGTH(nodes);
char *t = (char *)CHAR(STRING_ELT(target, 0));
int *status = NULL, *a = INTEGER(amat);

  if (isTRUE(debug))
    Rprintf("* (re)building cached information about node %s.\n", t);

  /* allocate and initialize the status vector. */
  status = alloc1dcont(length_nodes);

  /* iterate fo find the node position in the array.  */
  for (i = 0; i < length_nodes; i++)
    if (!strcmp(t, CHAR(STRING_ELT(nodes, i))))
      break;

  /* return the corresponding part of the bn structure. */
  return cache_node_structure(i, nodes, a, length_nodes, status, debuglevel);

}/*CACHE_PARTIAL_STRUCTURE*/

/* compute the cached values for a single node (C-friendly). */
SEXP c_cache_partial_structure(int target, SEXP nodes, SEXP amat, int *status, SEXP debug) {

int debuglevel = LOGICAL(debug)[0], length_nodes = LENGTH(nodes);
int *a = INTEGER(amat);

  /* allocate and initialize the status vector. */
  if (!(*status))
    status = alloc1dcont(length_nodes);

  /* return the corresponding part of the bn structure. */
  return cache_node_structure(target, nodes, a, length_nodes, status, debuglevel);

}/*C_CACHE_PARTIAL_STRUCTURE*/

/* backend to compute the cached values fro a single node. */
static SEXP cache_node_structure(int cur, SEXP nodes, int *amat, int nrows,
    int *status, int debuglevel) {

int i = 0, j = 0;
int num_parents = 0, num_children = 0, num_neighbours = 0, num_blanket = 0;
SEXP structure, names, mb, nbr, children, parents;

  if (debuglevel > 0)
    Rprintf("* node %s.\n", NODE(cur));

  for (i = 0; i < nrows; i++) {

    if (amat[CMC(cur, i, nrows)] == 1) {

      if (amat[CMC(i, cur, nrows)] == 0) {

        /* if a[i,j] = 1 and a[j,i] = 0, then i -> j. */
        if (debuglevel > 0)
          Rprintf("  > found child %s.\n", NODE(i));

        status[i] = CHILD;

        /* check whether this child has any other parent. */
        for (j = 0; j < nrows; j++) {

          if ((amat[CMC(j, i, nrows)] == 1) && (amat[CMC(i, j, nrows)] == 0)
                && (j != cur)) {

            /* don't mark a neighbour as in the markov blanket. */
            if (status[j] <= 1) {

              status[j] = BLANKET;

              if (debuglevel > 0)
                Rprintf("  > found node %s in markov blanket.\n", NODE(j));

            }/*THEN*/

          }/*THEN*/

        }/*FOR*/

      }/*THEN*/
      else {

        /* if a[i,j] = 1 and a[j,i] = 1, then i -- j. */
        if (debuglevel > 0)
          Rprintf("  > found neighbour %s.\n", NODE(i));

        status[i] = NEIGHBOUR;

      }/*ELSE*/

    }/*THEN*/
    else {

      if (amat[CMC(i, cur, nrows)] == 1) {

        /* if a[i,j] = 0 and a[j,i] = 1, then i <- j. */
        if (debuglevel > 0)
          Rprintf("  > found parent %s.\n", NODE(i));

        status[i] = PARENT;

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  /* count how may nodes fall in each category. */
  for (i = 0; i < nrows; i++) {

    switch(status[i]) {

      case CHILD:
        /* a child is also a neighbour and belongs into the markov blanket. */
        num_children++;
        num_neighbours++;
        num_blanket++;
        break;
      case PARENT:
        /* the same goes for a parent. */
        num_parents++;
        num_neighbours++;
        num_blanket++;
        break;
      case NEIGHBOUR:
        /* it's not known if this is parent or a children, but it's certainly a neighbour. */
        num_neighbours++;
        num_blanket++;
        break;
      case BLANKET:
        num_blanket++;
        break;
      default:
        /* this node is not even in the markov blanket. */
        break;

    }/*SWITCH*/

  }/*FOR*/

  if (debuglevel > 0)
    Rprintf("  > node %s has %d parent(s), %d child(ren), %d neighbour(s) and %d nodes in the markov blanket.\n",
      NODE(cur), num_parents, num_children, num_neighbours, num_blanket);

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
  PROTECT(mb = allocVector(STRSXP, num_blanket));
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

