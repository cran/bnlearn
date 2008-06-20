#include <R.h>
#include <Rinternals.h>

/*
 *  Coordinate system for an upper triangular matrix:
 *
 *  (row - 1) * ncols - row * (row - 1) / 2
 *
 *  the first element are the standard row major
 *  order coordinates; the second element is an
 *  adjustment to account for the missing lower
 *  half of the matrix.
 *
 */

/*
 * Beware: this function is based on the assumption that
 * each arc is unique in the arc list; otherwise the
 * counter may be wrong, leading to false negatives.
 */

#define COORDS(x,y) (x - 1) * n + y - x * (x - 1) / 2
#define ARC(i,col) INTEGER(arcs)[i + col * nrows]

SEXP is_dag(SEXP arcs, SEXP nnodes) {

int i = 0;
int nrows = LENGTH(arcs)/2;
int n = INTEGER(nnodes)[0];
short int checklist[COORDS(n, n)];
SEXP res;

  /* initialize the checklist. */
  memset(checklist, '\0', sizeof(short int) * (COORDS(n, n)));

  /* allocate the result. */
  PROTECT(res = allocVector(LGLSXP, 1));
  LOGICAL(res)[0] = TRUE;

  for (i = 0; i < nrows; i++) {

    /*  
     *  if row > column, reverse the order of the coordinates
     *  to fall into the upper half of the virtual adjacency matrix.
     */
    if (ARC(i, 0) > ARC(i, 1)) {

      if (checklist[COORDS(ARC(i, 1), ARC(i, 0)) - 1] == 0) {

        /* this arc is no present in checklist; update it. */
        checklist[COORDS(ARC(i, 1), ARC(i, 0)) - 1] = 1;

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

    }/*THEN*/
    else {

      if (checklist[COORDS(ARC(i, 0), ARC(i, 1)) - 1] == 0) {

        /* this arc is no present in checklist; update it. */
        checklist[COORDS(ARC(i, 0), ARC(i, 1)) - 1] = 1;

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

    }/*ELSE*/

  }/*FOR*/

  UNPROTECT(1);

  return res;

}/*IS_DAG*/
