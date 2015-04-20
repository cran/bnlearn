#include "include/rcore.h"
#include "include/matrix.h"
#include "include/graph.h"

SEXP shd(SEXP learned, SEXP golden, SEXP debug) {

int i = 0, j = 0, c1 = 0, c2 = 0, shd = 0, nnodes = 0;
int *lrn = NULL, *ref = NULL, debuglevel = isTRUE(debug);
SEXP temp, nodes, l, r;

  /* get the labels of the nodes. */
  temp = getListElement(learned, "nodes");
  nodes = getAttrib(temp, R_NamesSymbol);
  nnodes = length(nodes);

  /* get the arcs of the learned network. */
  temp = getListElement(learned, "arcs");
  /* build the adjacency matrix. */
  PROTECT(l = arcs2amat(temp, nodes));
  lrn = INTEGER(l);

  /* get the arcs of the learned network. */
  temp = getListElement(golden, "arcs");
  /* build the adjacency matrix. */
  PROTECT(r = arcs2amat(temp, nodes));
  ref = INTEGER(r);

  for (i = 0; i < nnodes; i++) {

    for (j = i + 1; j < nnodes; j++) {

      /* compute coordinates only once per iteration. */
      c1 = CMC(i, j, nnodes);
      c2 = CMC(j, i, nnodes);

      /* the two arcs are identical, nothing to do. */
      if ((lrn[c1] == ref[c1]) && (lrn[c2] == ref[c2]))
        continue;

      if (debuglevel > 0) {

        Rprintf("* arcs between %s and %s do not match.\n", NODE(i), NODE(j));

        if ((lrn[c1] == 1) && (lrn[c2] == 1))
          Rprintf("  > the learned network contains %s - %s.\n", NODE(i), NODE(j));
        else if ((lrn[c1] == 0) && (lrn[c2] == 0))
          Rprintf("  > the learned network contains no arc between %s and %s.\n", NODE(i), NODE(j));
        else if ((lrn[c1] == 1) && (lrn[c2] == 0))
          Rprintf("  > the learned network contains %s -> %s.\n", NODE(i), NODE(j));
        else if ((lrn[c1] == 0) && (lrn[c2] == 1))
          Rprintf("  > the learned network contains %s -> %s.\n", NODE(j), NODE(i));

        if ((ref[c1] == 1) && (ref[c2] == 1))
          Rprintf("  > the true network contains %s - %s.\n", NODE(i), NODE(j));
        else if ((ref[c1] == 0) && (ref[c2] == 0))
          Rprintf("  > the true network contains no arc between %s and %s.\n", NODE(i), NODE(j));
        else if ((ref[c1] == 1) && (ref[c2] == 0))
          Rprintf("  > the true network contains %s -> %s.\n", NODE(i), NODE(j));
        else if ((ref[c1] == 0) && (ref[c2] == 1))
          Rprintf("  > the true network contains %s -> %s.\n", NODE(j), NODE(i));

      }/*THEN*/

      /* increase the distance by one. */
      shd++;

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

  return ScalarInteger(shd);

}/*SHD*/

