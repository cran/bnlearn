#include "include/rcore.h"
#include "include/matrix.h"

/* Beware: this function is based on the assumption that each arc is unique
 * in the arc set; otherwise the counter may be wrong, leading to false
 * negatives. */

/* determine whether a graph is DAG or a PDAG/UG. */
SEXP is_dag(SEXP arcs, SEXP nodes) {

int i = 0, nrow = length(arcs)/2, n = LENGTH(nodes);
int *a = NULL;
short int *checklist = NULL;
SEXP try;

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  a = INTEGER(try);

  /* allocate and initialize the checklist. */
  checklist = Calloc1D(UPTRI_MATRIX(n), sizeof(short int));

  for (i = 0; i < nrow; i++) {

    if (checklist[UPTRI(a[CMC(i, 0, nrow)], a[CMC(i, 1, nrow)], n)] == 0) {

      /* this arc is not present in the checklist; add it. */
      checklist[UPTRI(a[CMC(i, 0, nrow)], a[CMC(i, 1, nrow)], n)] = 1;

    }/*THEN*/
    else {

      /* this arc or its opposite already present in the checklist; the graph
       * has at least an undirected arc, so return FALSE. */
      UNPROTECT(1);

      Free1D(checklist);

      return ScalarLogical(FALSE);

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  Free1D(checklist);

  return ScalarLogical(TRUE);

}/*IS_DAG*/

