#include "common.h"

/* return the complete orientation of a graph (the nodes argument gives
  * the node ordering). */
SEXP pdag2dag(SEXP arcs, SEXP nodes) {

int i = 0, j = 0, n = LENGTH(nodes);
int *a = NULL;
SEXP amat, res;

  /* build the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* scan the adjacency matrix. */
  for (i = 0; i < n; i++) {

    for (j = i + 1; j < n; j++) {

      /* if an arc is undirected, kill the orientation that violates the
       * specified node ordering (the one which is located in the lower
       * half of the matrix). */
      if ((a[CMC(i, j, n)] == 1) && (a[CMC(j, i, n)] == 1))
        a[CMC(j, i, n)] = 0;

    }/*FOR*/

  }/*FOR*/

  /* build the return value. */
  PROTECT(res = amat2arcs(amat, nodes));

  UNPROTECT(2);

  return res;

}/*PDAG2DAG*/

/* return the skeleton of a DAG/PDAG. */
SEXP dag2ug(SEXP arcs, SEXP nodes) {

int i = 0, j = 0, k = 0, coords = 0, narcs = LENGTH(arcs)/2, n = LENGTH(nodes);
int *a = NULL;
short int *checklist = NULL;
SEXP try, res, dimnames, colnames;

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  a = INTEGER(try);

  /* initialize the checklist. */
  checklist = allocstatus(UPTRI_MATRIX(n));

  /* add the arcs into the checklist. */
  for (i = 0; i < narcs; i++) {

      /* compute the index of the arc for the checklist array (which is the
       * upper-half of the adjacency matrix flattened into a 1-dimensional
       * array). */
      coords = UPTRI(a[CMC(i, 0, narcs)], a[CMC(i, 1, narcs)], n);

      /* count the number of arcs ignoring their orientation. */
      if (checklist[coords] == 0)
        k++;

      checklist[coords] = 1;

  }/*FOR*/

  /* unprotect try. */
  UNPROTECT(1);

  /* allocate the new arc set. */
  narcs = 2 * k;
  PROTECT(res = allocMatrix(STRSXP, narcs, 2));

  /* fill the arc set from the adjacency matrix. */
  for (i = 0, k = 0; i < n; i++) {

    for (j = i; j < n; j++) {

      coords = UPTRI(i + 1, j + 1, n);

      /* save both orientations for each arc. */
      if (checklist[coords] == 1) {

         SET_STRING_ELT(res, k, STRING_ELT(nodes, i));
         SET_STRING_ELT(res, k + 1 * narcs, STRING_ELT(nodes, j));
         k++;
         SET_STRING_ELT(res, k, STRING_ELT(nodes, j));
         SET_STRING_ELT(res, k + 1 * narcs, STRING_ELT(nodes, i));
         k++;

      }/*THEN*/

      /* no more arcs, get out of both loops. */
      if (k == narcs) goto end;

    }/*FOR*/

  }/*FOR*/

end:

  /* allocate colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  /* set the column names. */
  setAttrib(res, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);

  return res;

}/*DAG2UG*/

/* transform a DAG/PDAG/UDAG in a DAG with arc directions agreeing with the
 * node ordering specified by the "nodes" argmument. */
SEXP unique_arcs(SEXP arcs, SEXP nodes) {

int i = 0, j = 0, k = 0, nrows = 0, n = LENGTH(nodes);
short int *checklist = NULL;
int *coords = NULL;
SEXP result, dimnames, colnames, try;

  if (isNull(arcs)) {

    /* use NULL as a special jolly value which returns all possible arcs
     * given the specified node ordering. */
    nrows = n * (n - 1)/2;

    /* allocate the return value. */
    PROTECT(result = allocMatrix(STRSXP, nrows, 2));

    /* fill in the nodes' labels. */
    for (i = 0; i < n; i++) {

      for (j = i + 1; j < n; j++) {

        SET_STRING_ELT(result, CMC(k, 0, nrows), STRING_ELT(nodes, i));
        SET_STRING_ELT(result, CMC(k, 1, nrows), STRING_ELT(nodes, j));
        k++;

      }/*FOR*/

    }/*FOR*/

  }/*THEN*/
  else if (LENGTH(arcs) == 0) {

    /* the arc set is empty, nothing to do. */
    PROTECT(result = duplicate(arcs));

  }/*THEN*/
  else {

    /* there really is a non-empty arc set, process it. */
    nrows = LENGTH(arcs)/2;

    /* initialize the checklist. */
    checklist = allocstatus(UPTRI_MATRIX(n));
    
    /* match the node labels in the arc set. */
    PROTECT(try = match(nodes, arcs, 0));
    coords = INTEGER(try);

    /* indentify the arcs in the reduced adjacency matrix. */
    for (i = 0; i < nrows; i++)
      checklist[UPTRI(coords[CMC(i, 0, nrows)], coords[CMC(i, 1, nrows)], n)]++;

    UNPROTECT(1);

    /* count them; allocate and initialize the return value. */
    nrows = 0;

    for (i = 0; i < UPTRI_MATRIX(n); i++)
      if (checklist[i] > 0)
        nrows++;

    PROTECT(result = allocMatrix(STRSXP, nrows, 2));

    /* store the correct arcs in the return value. */
    for (i = 0; i < n; i ++) {

      for (j = i + 1; j < n; j++) {

        if (checklist[UPTRI(i + 1, j + 1, n)] > 0) {

          SET_STRING_ELT(result, CMC(k, 0, nrows), STRING_ELT(nodes, i));
          SET_STRING_ELT(result, CMC(k, 1, nrows), STRING_ELT(nodes, j));
          k++;

        }/*THEN*/

      }/*FOR*/

    }/*FOR*/

  }/*ELSE*/

  /* allocate, initialize and set the column names. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(result, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);

  return result;

}/*PDAG2DAG*/
