#include "include/rcore.h"
#include "include/graph.h"
#include "include/matrix.h"
#include "include/dataframe.h"
#include "include/globals.h"

/* adjusted arc counting for boot.strength(). */
SEXP bootstrap_strength_counters(SEXP prob, SEXP weight, SEXP arcs, SEXP nodes) {

int i = 0, j = 0, n = length(nodes), *a = NULL;
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

int i = 0, j = 0, k = 0, narcs = 0, nnodes = length(nodes);
double *p = NULL, *s = NULL, *d = NULL, tol = MACHINE_TOL;;
SEXP res, rownames, from, to, str, dir;

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
  setAttrib(res, R_ClassSymbol, mkString("data.frame"));

  /* allocate, initialize and set row names. */
  PROTECT(rownames = allocVector(INTSXP, narcs));
  for (i = 0; i < narcs; i++)
    INTEGER(rownames)[i] = i + 1;
  setAttrib(res, R_RowNamesSymbol, rownames);

  /* set column names. */
  setAttrib(res, R_NamesSymbol,
    mkStringVec(4, "from", "to", "strength", "direction"));

  /* attach the four columns. */
  SET_VECTOR_ELT(res, 0, from);
  SET_VECTOR_ELT(res, 1, to);
  SET_VECTOR_ELT(res, 2, str);
  SET_VECTOR_ELT(res, 3, dir);

  UNPROTECT(6);
  return res;

}/*BOOTSTRAP_ARC_COEFFICIENTS*/

/* reduce multiple boostrap strength R objects. */
SEXP bootstrap_reduce(SEXP x) {

int i = 0, j = 0, reps = length(x), nrows = 0;
double *str = NULL, *dir = NULL, *temp = NULL;
SEXP result, df, strength, direction;

  /* allocate return value. */
  PROTECT(result = allocVector(VECSXP, 4));

  /* extract the first data frame from the list. */
  df = VECTOR_ELT(x, 0);
  /* copy data frame column names. */
  setAttrib(result, R_NamesSymbol, getAttrib(df, R_NamesSymbol));
  /* copy the first two columns. */
  SET_VECTOR_ELT(result, 0, VECTOR_ELT(df, 0));
  SET_VECTOR_ELT(result, 1, VECTOR_ELT(df, 1));
  /* get the number of rows. */
  nrows = length(VECTOR_ELT(df, 0));
  /* allocate the remaining two columns. */
  PROTECT(strength = allocVector(REALSXP, nrows));
  str = REAL(strength);
  PROTECT(direction = allocVector(REALSXP, nrows));
  dir = REAL(direction);
  /* just copy over strength and direction. */
  memcpy(str, REAL(VECTOR_ELT(df, 2)), nrows * sizeof(double));
  memcpy(dir, REAL(VECTOR_ELT(df, 3)), nrows * sizeof(double));

  for (i = 1; i < reps; i++) {

    /* extract the data frame from the list. */
    df = VECTOR_ELT(x, i);
    /* accumulate strength. */
    temp = REAL(VECTOR_ELT(df, 2));
    for (j = 0; j < nrows; j++)
      str[j] += temp[j];
    /* accumulate direction. */
    temp = REAL(VECTOR_ELT(df, 3));
    for (j = 0; j < nrows; j++)
      dir[j] += temp[j];

  }/*FOR*/

  /* normalize dividing by the number of data frames. */
  for (j = 0; j < nrows; j++) {

    str[j] /= reps;
    dir[j] /= reps;

  }/*FOR*/

  /* set the last two columns. */
  SET_VECTOR_ELT(result, 2, strength);
  SET_VECTOR_ELT(result, 3, direction);
  /* make the return value a real data frame. */
  minimal_data_frame(result);

  UNPROTECT(3);

  return result;

}/*BOOTSTRAP_REDUCE*/

