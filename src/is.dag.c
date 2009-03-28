#include "common.h"

/*
 * Beware: these functions are based on the assumption
 * that each arc is unique in the arc set; otherwise
 * the counter may be wrong, leading to false negatives.
 */

#define ARC(i,col) INTEGER(arcs)[i + col * nrows]

SEXP is_dag(SEXP arcs, SEXP nnodes) {

int i = 0;
int nrows = LENGTH(arcs)/2;
int n = INT(nnodes);
short int *checklist;
SEXP res;

  /* allocate and initialize the checklist. */
  checklist = allocstatus(UPTRI(n, n, n));

  /* allocate the result. */
  PROTECT(res = allocVector(LGLSXP, 1));
  LOGICAL(res)[0] = TRUE;

  for (i = 0; i < nrows; i++) {

    if (checklist[UPTRI(ARC(i, 0), ARC(i, 1), n) - 1] == 0) {

      /* this arc is no present in checklist; update it. */
      checklist[UPTRI(ARC(i, 0), ARC(i, 1), n) - 1] = 1;

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

SEXP which_undirected(SEXP arcs) {

int i = 0;
int nrows = LENGTH(arcs)/2;
int n = LENGTH(getAttrib(arcs, R_LevelsSymbol));
short int *checklist;
SEXP res;

  /* initialize the checklist. */
  checklist = allocstatus(UPTRI(n, n, n));

  /* allocate the result. */
  PROTECT(res = allocVector(LGLSXP, nrows));

  for (i = 0; i < nrows; i++) {

    checklist[UPTRI(ARC(i, 0), ARC(i, 1), n) - 1]++;

  }/*FOR*/

  for (i = 0; i < nrows; i++) {

    if (checklist[UPTRI(ARC(i, 0), ARC(i, 1), n) - 1] == 1) {

      LOGICAL(res)[i] = FALSE;

    }/*THEN*/
    else {

      LOGICAL(res)[i] = TRUE;

    }/*ELSE*/

  }/*FOR*/

  UNPROTECT(1);

  return res;

}/*WHICH_UNDIRECTED*/
