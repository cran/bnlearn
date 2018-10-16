#include "include/rcore.h"
#include "include/sets.h"

/* initialize a one-dimensional contingency table. */
int fill_1d_table(int *xx, int **n, int llx, int num) {

int i = 0, ncomplete = 0;

  *n = Calloc1D(llx, sizeof(int));

  /* first fill the counts into the table... */
  for (i = 0; i < num; i++)
    if (xx[i] != NA_INTEGER)
      (*n)[xx[i] - 1]++;

  /* ... then add them up to count the number of complete observations. */
  for (i = 0; i < llx; i++)
    ncomplete += (*n)[i];

  return ncomplete;

}/*FILL_1D_TABLE*/

/* initialize a two-dimensional contingency table and the marginals. */
int fill_2d_table(int *xx, int *yy, int ***n, int **ni, int **nj, int llx,
    int lly, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  *n = (int **) Calloc2D(llx, lly, sizeof(int));
  *ni = (int *) Calloc1D(llx, sizeof(int));
  *nj = (int *) Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    if ((xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*n)[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

    (*ni)[i] += (*n)[i][j];
    (*nj)[j] += (*n)[i][j];

  }/*FOR*/

  /* compute the number of complete observations. */
  for (i = 0; i < llx; i++)
    ncomplete += (*ni)[i];

  return ncomplete;

}/*FILL_2D_TABLE*/

/* initialize a three-dimensional contingency table and the marginals. */
int fill_3d_table(int *xx, int *yy, int *zz, int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  *n = (int ***) Calloc3D(llz, llx, lly, sizeof(int));
  *ni = (int **) Calloc2D(llz, llx, sizeof(int));
  *nj = (int **) Calloc2D(llz, lly, sizeof(int));
  *nk = (int *) Calloc1D(llz, sizeof(int));

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    if ((zz[k] != NA_INTEGER) && (xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*n)[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        (*ni)[k][i] += (*n)[k][i][j];
        (*nj)[k][j] += (*n)[k][i][j];
        (*nk)[k] += (*n)[k][i][j];

      }/*FOR*/

  /* compute the number of complete observations. */
  for (k = 0; k < llz; k++)
    ncomplete += (*nk)[k];

  return ncomplete;

}/*FILL_3D_TABLE*/

/* minimal implementation of table(). */
SEXP minimal_table(SEXP dataframe, SEXP missing) {

int i = 0, nrow = length(VECTOR_ELT(dataframe, 0)), ncol = length(dataframe);
int *dd = NULL, *tt = NULL, **columns = NULL, *cfg = NULL;
double ncells = 1;
SEXP table, dims, dimnames, cur;

  /* prepare the dimensions. */
  PROTECT(dims = allocVector(INTSXP, ncol));
  dd = INTEGER(dims);
  PROTECT(dimnames = allocVector(VECSXP, ncol));
  setAttrib(dimnames, R_NamesSymbol, getAttrib(dataframe, R_NamesSymbol));

  /* dereference the data frame, extract the levels. */
  columns = (int **) Calloc1D(ncol, sizeof(int *));

  for (i = 0; i < ncol ; i++) {

    /* extract the column from the data frame... */
    cur = VECTOR_ELT(dataframe, i);
    /* ... dereference it... */
    columns[i] = INTEGER(cur);
    /* ... extract the number of levels... */
    dd[i] = NLEVELS(cur);
    /* ... and save them in the table dimensions. */
    SET_VECTOR_ELT(dimnames, i, getAttrib(cur, R_LevelsSymbol));

    ncells *= dd[i];

  }/*FOR*/

  if (ncells > INT_MAX) {

    Free1D(columns);

    UNPROTECT(2);

    error("attempting to create a table with more than INT_MAX cells.");

  }/*THEN*/

  /* allocate and dereference the table. */
  PROTECT(table = allocVector(INTSXP, ncells));
  tt = INTEGER(table);
  memset(tt, '\0', ncells * sizeof(int));

  /* prepare the configurations. */
  cfg = Calloc1D(nrow, sizeof(int));
  c_fast_config(columns, nrow, ncol, dd, cfg, NULL, 0);

  if (isTRUE(missing)) {

     for (i = 0; i < nrow; i++)
       if (cfg[i] != NA_INTEGER)
         tt[cfg[i]]++;

  }/*THEN*/
  else {

    for (i = 0; i < nrow; i++)
      tt[cfg[i]]++;

  }/*ELSE*/

  /* set the attributess for class and dimensions. */
  setAttrib(table, R_ClassSymbol, mkString("table"));
  setAttrib(table, R_DimSymbol, dims);
  setAttrib(table, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);

  Free1D(columns);
  Free1D(cfg);

  return table;

}/*MINIMAL_TABLE*/

