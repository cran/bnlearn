#include "include/rcore.h"
#include "include/graph.h"

/* convert a set of neighbourhoods into an arc set. */
SEXP nbr2arcs(SEXP nbr) {

int i = 0, j = 0, k = 0, narcs = 0;
int length_names = 0;
SEXP arcs, temp, names;

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

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    /* allocate an empty arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));

    UNPROTECT(1);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));

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

  UNPROTECT(1);

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

