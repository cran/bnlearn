#include "include/rcore.h"
#include "include/bn.h"

/* get root or leaf nodes of the graph. */
SEXP root_nodes(SEXP bn, SEXP leaves) {

short int *status = NULL;
int *get_leaves = INTEGER(leaves);
int i = 0, k = 0, counter = 0;
SEXP temp, temp2, nodes, node_data, labels, result;

  /* get to the nodes' data. */
  nodes = getListElement(bn, "nodes");
  /* this is for "bn.fit" objects. */
  if (isNull(nodes))
    nodes = bn;
  /* get the nodes' labels. */
  labels = getAttrib(nodes, R_NamesSymbol);
  /* allocate and initialize a status vector. */
  status = Calloc1D(length(nodes), sizeof(short int));

  for (i = 0; i < length(nodes); i++) {

    /* get the parents/children of this node. */
    node_data = VECTOR_ELT(nodes, i);

    if (*get_leaves == FALSE)
      temp = getListElement(node_data, "parents");
    else
      temp = getListElement(node_data, "children");

    /* this is not a root/leaf node, go on. */
    if (length(temp) != 0)
      continue;

    /* this takes care of dubious neighbours in "bn" objects. */
    temp = getListElement(node_data, "nbr");

    if (!isNull(temp)) {

      if (*get_leaves == FALSE)
        temp2 = getListElement(node_data, "children");
      else
        temp2 = getListElement(node_data, "parents");

      /* in partially directed graphs not all neighbours can be classified as
       * parents or children due to undirected arcs; return only nodes which
       * are not incident on any undirected arc. */
      if (length(temp) != length(temp2))
        continue;

    }/*THEN*/

    /* this is a root/leaf node, all right. */
    status[i] = 1;
    /* increase the counter. */
    counter++;

  }/*FOR*/

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(STRSXP, counter));

  for (i = 0; i < length(nodes); i++)
    if (status[i] == 1)
      SET_STRING_ELT(result, k++, STRING_ELT(labels, i));

  UNPROTECT(1);

  Free1D(status);

  return result;

}/*ROOT_NODES*/

/* build the arc set out of a "bn.fit" object. */
SEXP fit2arcs(SEXP bn) {

int i = 0, j = 0, k = 0, narcs = 0;
SEXP labels, node_data, children, result;

  /* get the nodes' labels. */
  labels = getAttrib(bn, R_NamesSymbol);

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

  UNPROTECT(1);

  return result;

}/*FIT2ARCS*/

#define MATCH_NODES(which, value) \
  temp = getListElement(node_data, which);     \
  PROTECT(try = match(labels, temp, 0));       \
  matched = INTEGER(try);                      \
  for (i = 0; i < length(try); i++) {          \
    if (status[matched[i] - 1] == 0) {         \
      status[matched[i] - 1] = value;          \
      counter++;                               \
     }/*THEN*/                                 \
  }/*FOR*/                                     \
  UNPROTECT(1);

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
  labels = getAttrib(bn, R_NamesSymbol);
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

  UNPROTECT(1);

  Free1D(status);

  return mb;

}/*FITTED_MB*/

/* return the size of the arc set. */
SEXP num_arcs(SEXP bn) {

int i = 0, is_fitted = 0, res = 0;
char *element = NULL;
SEXP nodes, node_data, temp;

  /* get to the nodes' data. */
  nodes = getListElement(bn, "nodes");
  /* check whether this is a "bn.fit" or "bn" object. */
  is_fitted = isNull(nodes);
  /* set the parameters for the object structure. */
  if (is_fitted) {

    nodes = bn;
    element = "parents";

  }/*THEN*/
  else {

   element = "nbr";

  }/*ELSE*/

  for (i = 0; i < length(nodes); i++) {

    /* get the parents/children of this node. */
    node_data = VECTOR_ELT(nodes, i);
    temp = getListElement(node_data, element);
    res += length(temp);

  }/*FOR*/

  /* summing up the neighbours counts each arc twice, dedeuplicate. */
  if (!is_fitted)
    res /= 2;

  return ScalarInteger(res);

}/*NUM_ARCS*/

