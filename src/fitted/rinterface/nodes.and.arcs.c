#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../minimal/common.h"
#include "../../minimal/strings.h"

/* build the arc set out of a "bn.fit" object. */
SEXP fit2arcs(SEXP bn) {

int i = 0, j = 0, k = 0, narcs = 0;
SEXP labels, node_data, children, result;

  /* get the nodes' labels. */
  PROTECT(labels = getAttrib(bn, R_NamesSymbol));

  /* first pass: count the number of arcs. */
  for (i = 0; i <  length(bn); i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* count its children. */
    narcs += length(getListElement(node_data, "children"));

  }/*FOR*/

  /* allocate the arc set. */
  PROTECT(result = allocMatrix(STRSXP, narcs, 2));
  /* set the column names. */
  setDimNames(result, R_NilValue, mkStringVec(2, "from", "to"));

  /* second pass: initialize the return value. */
  for (i = 0; i <  length(bn); i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* get its children. */
    children = getListElement(node_data, "children");

    for (j = 0; j < length(children); j++) {

      /* set the labels of the nodes incident on the arc. */
      SET_STRING_ELT(result, k, STRING_ELT(labels, i));
      SET_STRING_ELT(result, k + narcs, STRING_ELT(children, j));
      /* go to the next arc. */
      k++;

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

  return result;

}/*FIT2ARCS*/

#define MATCH_NODES(which, value) \
  do { \
  temp = getListElement(node_data, which); \
  PROTECT(try = match(labels, temp, 0)); \
  matched = INTEGER(try); \
  for (i = 0; i < length(try); i++) { \
    if (status[matched[i] - 1] == 0) { \
      status[matched[i] - 1] = value; \
      counter++; \
     }/*THEN*/ \
  }/*FOR*/ \
  UNPROTECT(1); \
  } while (0)

#define BLANKET		 1
#define NEIGHBOUR	 2
#define PARENT		 3
#define CHILD		 4
#define TARGET		 5

/* get the Markov blanket of a node from a fitted model. */
SEXP fitted_mb(SEXP bn, SEXP target) {

int i = 0, j = 0, target_node = 0, nnodes = 0, counter = 0;
int *matched = NULL;
short int *status = NULL;
SEXP mb, labels, try, temp, node_data;

  /* get the nodes' labels. */
  PROTECT(labels = getAttrib(bn, R_NamesSymbol));
  nnodes = length(labels);
  /* allocate and initialize a status vector. */
  status = Calloc1D(nnodes, sizeof(short int));

  /* match the label of the target node. */
  PROTECT(try = match(labels, target, 0));
  target_node = INT(try);
  UNPROTECT(1);

  /* mark the target node as such. */
  status[target_node - 1] = TARGET;

  /* match the parents and the children of the target node. */
  node_data = VECTOR_ELT(bn, target_node - 1);
  MATCH_NODES("parents", PARENT);
  MATCH_NODES("children", CHILD);

  /* now match the parents of each child. */
  for (j = 0; j < nnodes; j++) {

    /* this is not a child, go on. */
    if (status[j] != CHILD)
      continue;

    node_data = VECTOR_ELT(bn, j);
    MATCH_NODES("parents", BLANKET);

  }/*FOR*/

  /* a node is not considered part of its own Markov blanket. */
  status[target_node - 1] = 0;

  /* allocate and initialize the result. */
  PROTECT(mb = allocVector(STRSXP, counter));
  for (i = 0, j = 0; i < nnodes; i++)
    if (status[i] != 0)
      SET_STRING_ELT(mb, j++, STRING_ELT(labels, i));

  UNPROTECT(2);

  Free1D(status);

  return mb;

}/*FITTED_MB*/

