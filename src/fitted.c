#include "common.h"

#define ROOTNODES 0

/* get root or leaf nodes of the graph. */
SEXP root_nodes(SEXP bn, SEXP check) {

short int *status = NULL;
int *to_check = INTEGER(check);
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
  status = allocstatus(LENGTH(nodes));

  for (i = 0; i < LENGTH(nodes); i++) {

    /* get the parents/children of this node. */
    node_data = VECTOR_ELT(nodes, i);

    if (*to_check == ROOTNODES)
      temp = getListElement(node_data, "parents");
    else
      temp = getListElement(node_data, "children");

    /* this is not a root/leaf node, go on. */
    if (LENGTH(temp) != 0)
      continue;

    /* this takes care of dubious neighbours in "bn" objects. */
    temp = getListElement(node_data, "nbr");

    if (!isNull(temp)) {

      if (*to_check == ROOTNODES)
        temp2 = getListElement(node_data, "children");
      else
        temp2 = getListElement(node_data, "parents");

      /* in partially directed graphs not all neighbours can be classified as
       * parents or children due to undirected arcs; return only nodes which
       * are not incident on any undirected arc. */
      if (LENGTH(temp) != LENGTH(temp2))
        continue;

    }/*THEN*/

    /* this is a root/leaf node, all right. */
    status[i] = 1;
    /* increase the counter. */
    counter++;

  }/*FOR*/

  /* allocate and initialize the result. */
  PROTECT(result = allocVector(STRSXP, counter));

  for (i = 0; i < LENGTH(nodes); i++)
    if (status[i] == 1)
      SET_STRING_ELT(result, k++, STRING_ELT(labels, i));

  UNPROTECT(1);

  return result;

}/*ROOT_NODES*/

/* build the arc set out of a "bn.fit" object. */
SEXP fit2arcs(SEXP bn) {

int i = 0, j = 0, k = 0, narcs = 0;
SEXP labels, node_data, children, result, colnames, dimnames;

  /* get the nodes' labels. */
  labels = getAttrib(bn, R_NamesSymbol);

  /* first pass: count the number of arcs. */
  for (i = 0; i <  LENGTH(bn); i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* count its children. */
    narcs += LENGTH(getListElement(node_data, "children"));

  }/*FOR*/

  /* allocate the column names. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* allocate the arc set. */
  PROTECT(result = allocMatrix(STRSXP, narcs, 2));
  /* set the column names. */
  setAttrib(result, R_DimNamesSymbol, dimnames);

  /* second pass: initialize the return value. */
  for (i = 0; i <  LENGTH(bn); i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* get its children. */
    children = getListElement(node_data, "children");

    for (j = 0; j < LENGTH(children); j++) {

      /* set the labels of the nodes incident on the arc. */
      SET_STRING_ELT(result, k, STRING_ELT(labels, i));
      SET_STRING_ELT(result, k + narcs, STRING_ELT(children, j));
      /* go to the next arc. */
      k++;

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(3);

  return result;

}/*FIT2ARCS*/

/* compute the number of parameters of the model. */
SEXP fitted_nparams(SEXP bn, SEXP debug) {

int i = 0, j = 0, node_params = 0, nnodes = LENGTH(bn);
int *res = NULL, *debuglevel = LOGICAL(debug);
SEXP result, nodes = R_NilValue, node_data, temp;

  /* allocate, dereference and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, 1));
  res = INTEGER(result);
  res[0] = 0;

  if (*debuglevel > 0)
    nodes = getAttrib(bn, R_NamesSymbol);

  for (i = 0; i < nnodes; i++) {

    /* get the node's data. */
    node_data = VECTOR_ELT(bn, i);
    /* get its probability distribution (if discrete). */
    temp = getListElement(node_data, "prob");

    if (!isNull(temp)) {

      /* reset the parameters' counter for this node. */
      node_params = 1;
      /* get the dimensions of the conditional probability table. */
      temp = getAttrib(temp, R_DimSymbol);
      /* compute the number of parameters. */
      for (j = 1; j < LENGTH(temp); j++)
        node_params *= INTEGER(temp)[j];

      node_params *= INTEGER(temp)[0] - 1;

    }/*THEN*/
    else {

      /* this is a continuous node, so it's a lot easier. */
      node_params = LENGTH(getListElement(node_data, "coefficients"));

    }/*ELSE*/

    if (*debuglevel > 0)
      Rprintf("* node %s has %d parameter(s).\n", NODE(i), node_params);

    res[0] += node_params;

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*FITTED_NPARAMS*/

#define MATCH_NODES(which, value) \
  temp = getListElement(node_data, which);     \
                                               \
  PROTECT(try = match(labels, temp, 0));       \
  matched = INTEGER(try);                      \
                                               \
  for (i = 0; i < LENGTH(try); i++)            \
    if (status[matched[i] - 1] == 0) {         \
                                               \
      status[matched[i] - 1] = value;          \
      counter++;                               \
                                               \
     }                                         \
                                               \
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
  nnodes = LENGTH(labels);
  /* allocate and initialize a status vector. */
  status = allocstatus(nnodes);

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

  return mb;

}/*FITTED_MB*/
