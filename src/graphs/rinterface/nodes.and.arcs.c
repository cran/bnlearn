#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../minimal/common.h"

/* get root or leaf nodes of the graph. */
SEXP root_nodes(SEXP bn, SEXP leaves) {

short int *status = NULL;
int i = 0, k = 0, counter = 0, get_leaves = isTRUE(leaves);
SEXP temp, temp2, nodes, node_data, labels, result;

  /* get to the nodes' data. */
  nodes = getListElement(bn, "nodes");
  /* this is for "bn.fit" objects. */
  if (isNull(nodes))
    nodes = bn;
  /* get the nodes' labels. */
  PROTECT(labels = getAttrib(nodes, R_NamesSymbol));
  /* allocate and initialize a status vector. */
  status = Calloc1D(length(nodes), sizeof(short int));

  for (i = 0; i < length(nodes); i++) {

    /* get the parents/children of this node. */
    node_data = VECTOR_ELT(nodes, i);

    if (get_leaves == FALSE)
      temp = getListElement(node_data, "parents");
    else
      temp = getListElement(node_data, "children");

    /* this is not a root/leaf node, go on. */
    if (length(temp) != 0)
      continue;

    /* this takes care of dubious neighbours in "bn" objects. */
    temp = getListElement(node_data, "nbr");

    if (!isNull(temp)) {

      if (get_leaves == FALSE)
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

  UNPROTECT(2);

  Free1D(status);

  return result;

}/*ROOT_NODES*/

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
