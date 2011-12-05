#include "common.h"

/* adjusted arc counting for boot.strength(). */
SEXP bootstrap_strength_counters(SEXP prob, SEXP weight, SEXP arcs, SEXP nodes) {

int i = 0, j = 0, n = LENGTH(nodes), *a = NULL;
double *p = NULL, *w = NULL;
SEXP amat;

  /* build the adjacency matrix for the current network. */
  PROTECT(amat = arcs2amat(arcs, nodes));

  /* map the contents of the SEXPs for easy access.  */
  a = INTEGER(amat);
  p = REAL(prob);
  w = REAL(weight);

  for (i = 0; i < n; i++) {

    for (j = 0; j < n; j++) {

      /* increase the counter of 1/2 for an undirected arc (the other half
       * is added to the symmetric element in the matrix) or of 1 for a
       * direcxted arc. */
      if (a[CMC(i, j, n)] == 1) {

        if (a[CMC(j, i, n)] == 1)
          p[CMC(i, j, n)] += 0.5 * (*w);
        else
          p[CMC(i, j, n)] += 1 * (*w);

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);
  return prob;

}/*BOOTSTRAP_STRENGTH*/

/* arc strength (confidence) and direction coefficients. */
SEXP bootstrap_arc_coefficients(SEXP prob, SEXP nodes) {

int i = 0, j = 0, k = 0, narcs = 0, nnodes = LENGTH(nodes);
double *p = NULL, *s = NULL, *d = NULL, tol = MACHINE_TOL;;
SEXP res, class, rownames, colnames, from, to, str, dir;

  /* compute the dimension of the arcs set. */
  narcs = nnodes * (nnodes - 1);

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
      /* sanitize out-of-boundary values arising from floating point errors. */
      s[k] = (s[k] < tol) ? 0 : s[k];
      s[k] = (s[k] > 1 - tol) ? 1 : s[k];
      d[k] = (d[k] < tol) ? 0 : d[k];
      d[k] = (d[k] > 1 - tol) ? 1 : d[k];

      /* increment the arc counter. */
      k++;

    }/*FOR*/

  }/*FOR*/

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

