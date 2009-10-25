#include "common.h"

/* adjusted arc counting for boot.strength(). */
SEXP bootstrap_strength_counters(SEXP prob, SEXP arcs, SEXP nodes) {

int i = 0, j = 0, n = LENGTH(nodes), *a = NULL;
double *p;
SEXP amat;

  /* build the adjacency matrix for the current network. */
  PROTECT(amat = arcs2amat(arcs, nodes));

  /* map the contents of the SEXPs for easy access.  */
  a = INTEGER(amat);
  p = REAL(prob);

  for (i = 0; i < n; i++) {

    for (j = 0; j < n; j++) {

      /* increase the counter of 1/2 for an undirected arc (the other half
       * is added to the symmetric element in the matrix) or of 1 for a
       * direcxted arc. */
      if (a[CMC(i, j, n)] == 1) {

        if (a[CMC(j, i, n)] == 1)
          p[CMC(i, j, n)] += 0.5;
        else
          p[CMC(i, j, n)] += 1;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);
  return prob;

}/*BOOTSTRAP_STRENGTH*/

/* arc strength and direction confidence coefficients. */
SEXP bootstrap_arc_coefficients(SEXP prob, SEXP arcs, SEXP nodes) {

int i = 0, j = 0, k = 0, *coords = NULL, narcs = 0, nnodes = LENGTH(nodes);
double *p = NULL, *s = NULL, *d = NULL;
SEXP res, class, rownames, colnames, from, to, str, dir, try;

  /* get the dimension of the arcs set; if none is specified get them all. */
  narcs = (isNull(arcs) ? nnodes * nnodes - nnodes : LENGTH(arcs) / 2);

  /* allocate and initialize the various columns. */
  PROTECT(from = allocVector(STRSXP, narcs));
  PROTECT(to = allocVector(STRSXP, narcs));
  PROTECT(str = allocVector(REALSXP, narcs));
  PROTECT(dir = allocVector(REALSXP, narcs));

  /* dereference the probability matrix and the coefficients once and
   * for all. */
  p = REAL(prob);
  s = REAL(str);
  d = REAL(dir);

  if (isNull(arcs)) {

    /* fill in the coefficients. */
    for (i = 0, k = 0; i < nnodes; i++) {

      for (j = 0; j < nnodes; j++) {

        /* "from" must differ from "to". */
        if (i == j)
          continue;

        /* set the labels of the incident nodes. */
        SET_STRING_ELT(from, k, STRING_ELT(nodes, i));
        SET_STRING_ELT(to, k, STRING_ELT(nodes, j));
        /* compute arc strength and direction confidence. */
        s[k] = p[CMC(i, j, nnodes)] + p[CMC(j, i, nnodes)];
        d[k] = (s[k] == 0 ? 0 : p[CMC(i, j, nnodes)] / s[k]);

        /* increment the arc counter. */
        k++;

      }/*FOR*/

    }/*FOR*/

  }/*THEN*/
  else {

    /* match the node labels in the arc set. */
    PROTECT(try = match(nodes, arcs, 0));
    coords = INTEGER(try);

    for (k = 0; k < narcs; k++) {

      /* set the labels in the data frame. */
      SET_STRING_ELT(from, k, STRING_ELT(arcs, k));
      SET_STRING_ELT(to, k, STRING_ELT(arcs, k + narcs));

      /* match the positions of the nodes' labels and compute arc strength
       * and direction confidence of the arc. */
      s[k] = p[CMC(coords[k] - 1, coords[k + narcs] - 1, nnodes)] +
             p[CMC(coords[k + narcs] - 1, coords[k] - 1, nnodes)];
      d[k] = (s[k] == 0 ? 0 : p[CMC(coords[k], coords[k + narcs], nnodes)] / s[k]);

    }/*FOR*/

    UNPROTECT(1);

  }/*ELSE*/

  /* allocate and initialize the return value. */
  PROTECT(res = allocVector(VECSXP, 4));

  /* allocate, initialize and set the class name. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("data.frame"));
  setAttrib(res, R_ClassSymbol, class);

  /* allocate, initialize and set row names. */
  PROTECT(rownames = allocVector(INTSXP, narcs));
  for (i = 0; i < narcs; i++)
    INTEGER(rownames)[i] = i + 1;
  setAttrib(res, R_RowNamesSymbol, rownames);

  /* allocate, initialize and set column names. */
  PROTECT(colnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_STRING_ELT(colnames, 2, mkChar("strength"));
  SET_STRING_ELT(colnames, 3, mkChar("direction"));
  setAttrib(res, R_NamesSymbol, colnames);

  /* attach the four columns. */
  SET_VECTOR_ELT(res, 0, from);
  SET_VECTOR_ELT(res, 1, to);
  SET_VECTOR_ELT(res, 2, str);
  SET_VECTOR_ELT(res, 3, dir);

  UNPROTECT(8);
  return res;

}/*BOOTSTRAP_ARC_COEFFICIENTS*/

