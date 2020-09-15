#include "include/rcore.h"
#include "include/sets.h"
#include "include/data.structures.h"

/* create a one-dimensional contingency table. */
counts1d new_1d_table(int llx) {

counts1d table = { 0 };

  table.llx = llx;
  table.n = Calloc1D(llx, sizeof(int));

  return table;

}/*NEW_1D_TABLE*/

/* initialize a one-dimensional contingency table. */
void fill_1d_table(int *xx, counts1d *table, int num) {

int i = 0, ncomplete = 0;

  /* first fill the counts into the table... */
  for (i = 0; i < num; i++)
    if (xx[i] != NA_INTEGER)
      (*table).n[xx[i] - 1]++;

  /* ... then add them up to count the number of complete observations. */
  for (i = 0; i < (*table).llx; i++)
    ncomplete += (*table).n[i];
  (*table).nobs = ncomplete;

}/*FILL_1D_TABLE*/

/* print a one-dimensional contingency table. */
void print_1d_table(counts1d table) {

  Rprintf("1-dimensional contingency table (%d cells)\n", table.llx);

  for (int i = 0; i < table.llx; i++)
    Rprintf("%d ", table.n[i]);
  Rprintf("\n");

}/*PRINT_1D_TABLE*/

/* free a one-dimensional contingency table. */
void Free1DTAB(counts1d table) {

  Free1D(table.n);

}/*FREE1DTAB*/

/* create a two-dimensional contingency table. */
counts2d new_2d_table(int llx, int lly, bool margins) {

counts2d table = { 0 };

  table.llx = llx;
  table.lly = lly;
  table.n = (int **)Calloc2D(llx, lly, sizeof(int));

  if (margins) {

    table.ni = Calloc1D(llx, sizeof(int));
    table.nj = Calloc1D(lly, sizeof(int));

  }/*THEN*/

  return table;

}/*NEW_2D_TABLE*/

/* initialize a two-dimensional contingency table including the marginals. */
void fill_2d_table(int *xx, int *yy, counts2d *table, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    if ((xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*table).n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals if they have been allocated memory for. */
  if ((*table).ni && (*table).nj) {

    for (i = 0; i < (*table).llx; i++)
      for (j = 0; j < (*table).lly; j++) {

      (*table).ni[i] += (*table).n[i][j];
      (*table).nj[j] += (*table).n[i][j];

    }/*FOR*/

    /* compute the number of complete observations. */
    for (i = 0; i < (*table).llx; i++)
      ncomplete += (*table).ni[i];
    (*table).nobs = ncomplete;

  }/*THEN*/
  else {

    /* compute the number of complete observations. */
    for (i = 0; i < (*table).llx; i++)
      for (j = 0; j < (*table).lly; j++)
        ncomplete += (*table).n[i][j];
    (*table).nobs = ncomplete;

  }/*ELSE*/


}/*FILL_2D_TABLE*/

/* print a two-dimensional contingency table. */
void print_2d_table(counts2d table) {

  Rprintf("2-dimensional contingency table (%d x %d cells)\n",
    table.llx, table.lly);
  for (int i = 0; i < table.llx; i++) {

    for (int j = 0; j < table.lly; j++)
      Rprintf("%d ", table.n[i][j]);
    Rprintf("\n");

  }/*FOR*/

}/*PRINT_2D_TABLE*/

/* free a two-dimensional contingency table. */
void Free2DTAB(counts2d table) {

    Free2D(table.n, table.llx);
    Free1D(table.ni);
    Free1D(table.nj);

}/*FREE2DTAB*/

/* create a three-dimensional contingency table. */
counts3d new_3d_table(int llx, int lly, int llz) {

counts3d table = { 0 };

  table.llx = llx;
  table.lly = lly;
  table.llz = llz;
  table.n = (int ***) Calloc3D(llz, llx, lly, sizeof(int));
  table.ni = (int **) Calloc2D(llz, llx, sizeof(int));
  table.nj = (int **) Calloc2D(llz, lly, sizeof(int));
  table.nk = (int *) Calloc1D(llz, sizeof(int));

  return table;

}/*NEW_3D_TABLE*/

/* initialize a three-dimensional contingency table and the marginals. */
void fill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    if ((zz[k] != NA_INTEGER) && (xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*table).n[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < (*table).llx; i++)
    for (j = 0; j < (*table).lly; j++)
      for (k = 0; k < (*table).llz; k++) {

        (*table).ni[k][i] += (*table).n[k][i][j];
        (*table).nj[k][j] += (*table).n[k][i][j];
        (*table).nk[k] += (*table).n[k][i][j];

      }/*FOR*/

  /* compute the number of complete observations. */
  for (k = 0; k < (*table).llz; k++)
    ncomplete += (*table).nk[k];
  (*table).nobs = ncomplete;

}/*FILL_3D_TABLE*/

/* print a three-dimensional contingency table. */
void print_3d_table(counts3d table) {

  Rprintf("3-dimensional contingency table (%d x %d x %d cells)\n",
    table.llx, table.lly, table.llz);

  for (int k = 0; k < table.llz; k++) {

    Rprintf("[slice %d]", k);

    for (int i = 0; i < table.llx; i++) {

      for (int j = 0; j < table.lly; j++)
        Rprintf("%d ", table.n[k][i][j]);
      Rprintf("\n");

    }/*FOR*/

  }/*FOR*/

}/*PRINT_3D_TABLE*/

/* free a three-dimensional contingency table. */
void Free3DTAB(counts3d table) {

  Free3D(table.n, table.llz, table.llx);
  Free2D(table.ni, table.llz);
  Free2D(table.nj, table.llz);
  Free1D(table.nk);

}/*FREE3DTAB*/

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

