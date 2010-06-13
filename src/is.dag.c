#include "common.h"

/*
 *  Beware: this function is based on the assumption that each arc is unique
 *  in the arc set; otherwise the counter may be wrong, leading to false
 *  negatives.
 */

/* determine whether a graph is DAG or a PDAG/UG. */
SEXP is_dag(SEXP arcs, SEXP nnodes) {

int i = 0, nrows = LENGTH(arcs)/2, n = INT(nnodes);
int *a = INTEGER(arcs);
short int *checklist = NULL;
SEXP res;

  /* allocate and initialize the checklist. */
  checklist = allocstatus(UPTRI_MATRIX(n));

  /* allocate the result. */
  PROTECT(res = allocVector(LGLSXP, 1));
  LOGICAL(res)[0] = TRUE;

  for (i = 0; i < nrows; i++) {

    if (checklist[UPTRI(a[CMC(i, 0, nrows)], a[CMC(i, 1, nrows)], n)] == 0) {

      /* this arc is not present in the checklist; add it. */
      checklist[UPTRI(a[CMC(i, 0, nrows)], a[CMC(i, 1, nrows)], n)] = 1;

    }/*THEN*/
    else {

      /*  this arc or its opposite already present in the
       *  checklist; the graph has at least an undirected
       *  arc, so return FALSE.
       */
      LOGICAL(res)[0] = FALSE;
      UNPROTECT(1);
      return res;

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  return res;

}/*IS_DAG*/

