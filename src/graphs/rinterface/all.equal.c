#include "../../include/rcore.h"
#include "../../core/sort.h"
#include "../../include/graph.h"
#include "../../minimal/common.h"

/* backend for the all.equal() function for bn objects. */
SEXP all_equal_bn(SEXP target, SEXP current) {

int nnodes = 0,  narcs = 0;
int *t = NULL, *c = NULL;
SEXP tnodes, cnodes, cmatch, tarcs, carcs, thash, chash;

  /* get the node set of each network. */
  PROTECT(tnodes = getAttrib(getListElement(target, "nodes"), R_NamesSymbol));
  PROTECT(cnodes = getAttrib(getListElement(current, "nodes"), R_NamesSymbol));

  /* first check: node sets must have the same size. */
  if (length(tnodes) != length(cnodes)) {

    UNPROTECT(2);

    return mkString("Different number of nodes");

  }/*THEN*/

  /* store for future use. */
  nnodes = length(tnodes);

  /* second check: node sets must contain the same node labels.  */
  PROTECT(cmatch = match(tnodes, cnodes, 0));
  c = INTEGER(cmatch);
  /* sorting takes care of different node orderings. */
  i_sort(c, NULL, nnodes);

  /* check that every node in the first network is also present in the
   * second one; this is enough because the node sets have the same size
   * and the nodes in each set are guaranteed to be unique. */
  for (int i = 0; i < nnodes; i++) {

    if (c[i] != i + 1) {

      UNPROTECT(3);

      return mkString("Different node sets");

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  /* get the node set of each network. */
  tarcs = getListElement(target, "arcs");
  carcs = getListElement(current, "arcs");

  /* third check: arc sets must have the same size. */
  if (length(tarcs) != length(carcs)) {

    UNPROTECT(2);

    return mkString("Different number of directed/undirected arcs");

  }/*THEN*/

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
    i_sort(t, NULL, narcs);
    i_sort(c, NULL, narcs);

    /* compare the integer vectors as generic memory areas. */
    if (memcmp(t, c, narcs * sizeof(int))) {

      UNPROTECT(4);

      return mkString("Different arc sets");

    }/*THEN*/

    UNPROTECT(2);

  }/*THEN*/

  UNPROTECT(2);

  /* all checks completed successfully, returning TRUE. */
  return ScalarLogical(TRUE);

}/*ALL_EQUAL_BN*/

