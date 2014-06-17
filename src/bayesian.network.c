#include "common.h"

/* get the number of parameters of the whole network (discrete case). */
SEXP nparams_dnet(SEXP graph, SEXP data, SEXP real, SEXP debug) {

int i = 0, j = 0, nnodes = 0;
int *index = NULL, r = isTRUE(real), debuglevel = isTRUE(debug);
double node_params = 0, all_params = 0;
double *nlevels = NULL;
SEXP nodes, node_data, parents, try;

  /* get nodes' number and data. */
  node_data = getListElement(graph, "nodes");
  nodes = getAttrib(node_data, R_NamesSymbol);
  nnodes = length(node_data);

  /* get the level count for each node. */
  nlevels = alloc1dreal(nnodes);
  for (i = 0; i < nnodes; i++)
    nlevels[i] = NLEVELS2(data, i);

  /* for each node... */
  for (i = 0; i < nnodes; i++) {

    /* reset the parameter counter. */
    node_params = 1;

    /* match the parents of the node. */
    parents = getListElement(VECTOR_ELT(node_data, i), "parents");
    PROTECT(try = match(nodes, parents, 0));
    index = INTEGER(try);

    /* compute the number of configurations. */
    for (j = 0; j < length(try); j++)
      node_params *= nlevels[index[j] - 1];

    UNPROTECT(1);

    /* multiply by the number of free parameters. */
    if (r > 0)
      node_params *= nlevels[i] - 1;
    else
      node_params *= nlevels[i];

    if (debuglevel > 0)
      Rprintf("* node %s has %.0lf parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

  }/*FOR*/

  return ScalarReal(all_params);

}/*NPARAMS_DNET*/

/* get the number of parameters of a single node (discrete case). */
SEXP nparams_dnode(SEXP graph, SEXP node, SEXP data, SEXP real) {

int i = 0, j = 0, length_nodes = 0, r = isTRUE(real);
double nlevels = 1;
SEXP temp, names;

  /* get the entry for the parents of the node.*/
  temp = getListElement(graph, "nodes");
  temp = getListElement(temp, (char *)CHAR(STRING_ELT(node, 0)));
  temp = getListElement(temp, "parents");

  /* get the column names from the data set and the length of the
       relevant vectors. */
  names = getAttrib(data, R_NamesSymbol);
  length_nodes = length(temp);

  /* sum (multiply, actually) up the levels. */
  for (i = 0; i < length(names); i++) {

    for (j = 0; j < length_nodes; j++) {

      /* this is a parent. */
      if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(temp, j)))) {

        nlevels *= NLEVELS2(data, i);

      }/*THEN*/

    }/*FOR*/

    /* this is the node. */
    if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(node, 0)))) {

      nlevels *= NLEVELS2(data, i) - 1 * r;

    }/*THEN*/

  }/*FOR*/

  return ScalarReal(nlevels);

}/*NPARAMS_DNODE*/

/* get the number of parameters of the whole network (continuous case). */
SEXP nparams_gnet(SEXP graph, SEXP debug) {

int i = 0, node_params = 0, debuglevel = isTRUE(debug);
double all_params = 0;
SEXP nodes = R_NilValue, temp = getListElement(graph, "nodes");

  if (debuglevel > 0)
    nodes = getAttrib(temp, R_NamesSymbol);

  /* add one parameter for each regressor, which means one for each
   * parent for each node plus the intercept. */
  for (i = 0; i < length(temp); i++) {

    node_params = length(getListElement(VECTOR_ELT(temp, i), "parents")) + 1;

    if (debuglevel > 0)
      Rprintf("* node %s has %d parameter(s).\n", NODE(i), node_params);

    /* update the return value. */
    all_params += node_params;

  }/*FOR*/

  return ScalarReal(all_params);

}/*NPARAMS_GNET*/

/* get the number of parameters of a single node (continuous case). */
SEXP nparams_gnode(SEXP graph, SEXP node) {

char *name = (char *)CHAR(STRING_ELT(node, 0));
SEXP temp;

  temp = getListElement(graph, "nodes");
  temp = getListElement(temp, name);
  temp = getListElement(temp, "parents");

  return ScalarReal(length(temp) + 1);

}/*NPARAMS_GNODE*/

/* convert a set of neighbourhoods into an arc set. */
SEXP nbr2arcs(SEXP nbr) {

int i = 0, j = 0, k = 0, narcs = 0;
int length_names = 0;
SEXP arcs, dimnames, temp, names;

  /* get the names of the nodes. */
  names = getAttrib(nbr, R_NamesSymbol);
  length_names = length(names);

  /* scan the structure to determine the number of arcs.  */
  for (i = 0; i < length_names; i++) {

    /* get the entry for the neighbours of the node.*/
    temp = getListElement(nbr, (char *)CHAR(STRING_ELT(names, i)));
    temp = getListElement(temp, "nbr");

    narcs += length(temp);

  }/*FOR*/

  /* allocate colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, mkStringVec(2, "from", "to"));

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    /* allocate an empty arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

  }/*ELSE*/

  /* rescan the structure to build the arc set. */
  for (i = 0; i < length_names; i++) {

    /* get the entry for the neighbours of the node.*/
    temp = getListElement(nbr, (char *)CHAR(STRING_ELT(names, i)));
    temp = getListElement(temp, "nbr");

    for (j = 0; j < length(temp); j++) {

      SET_STRING_ELT(arcs, k, STRING_ELT(names, i));
      SET_STRING_ELT(arcs, k + 1 * narcs , STRING_ELT(temp, j));
      k++;

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

  return arcs;

}/*NBR2ARCS*/

/* backend for the all.equal() function for bn objects. */
SEXP all_equal(SEXP target, SEXP current) {

int nnodes = 0,  narcs = 0;
int *t = NULL, *c = NULL;
SEXP tnodes, cnodes, cmatch, tarcs, carcs, thash, chash;

  /* get the node set of each network. */
  tnodes = getAttrib(getListElement(target, "nodes"), R_NamesSymbol);
  cnodes = getAttrib(getListElement(current, "nodes"), R_NamesSymbol);

  /* first check: node sets must have the same size. */
  if (length(tnodes) != length(cnodes))
    return mkString("Different number of nodes");

  /* store for future use. */
  nnodes = length(tnodes);

  /* second check: node sets must contain the same node labels.  */
  PROTECT(cmatch = match(tnodes, cnodes, 0));
  c = INTEGER(cmatch);
  /* sorting takes care of different node orderings. */
  R_isort(c, nnodes);

  /* check that every node in the first network is also present in the
   * second one; this is enough because the node sets have the same size
   * and the nodes in each set are guaranteed to be unique. */
  for (int i = 0; i < nnodes; i++) {

    if (c[i] != i + 1) {

      UNPROTECT(1);

      return mkString("Different node sets");

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  /* get the node set of each network. */
  tarcs = getListElement(target, "arcs");
  carcs = getListElement(current, "arcs");

  /* third check: arc sets must have the same size. */
  if (length(tarcs) != length(carcs))
    return mkString("Different number of directed/undirected arcs");

  /* store for future use. */
  narcs = length(tarcs)/2;

  /* fourth check:  arcs sets must contain the same arcs. */
  if (narcs > 0) {

    /* compute the numeric hashes of both arc sets (against the same
     * node set to make comparisons meaningful) and sort them. */
    PROTECT(thash = arc_hash(tarcs, tnodes, FALSE, TRUE));
    PROTECT(chash = arc_hash(carcs, tnodes, FALSE, TRUE));
    /* dereference the resulting integer vectors. */
    t = INTEGER(thash);
    c = INTEGER(chash);
    /* sorting takes care of different arc orderings. */
    R_isort(t, narcs);
    R_isort(c, narcs);

    /* compare the integer vectors as generic memory areas. */
    if (memcmp(t, c, narcs * sizeof(int))) {

      UNPROTECT(2);

      return mkString("Different arc sets");

    }/*THEN*/

    UNPROTECT(2);

  }/*THEN*/

  /* all checks completed successfully, returning TRUE. */
  return ScalarLogical(TRUE);

}/*ALL_EQUAL*/

