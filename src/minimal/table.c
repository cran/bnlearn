#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/sets.h"
#include "../math/linear.algebra.h"
#include "common.h"

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

  /* after this it is safe to cast ncels to an integer. */
  if (ncells > INT_MAX) {

    Free1D(columns);
    UNPROTECT(2);
    error("attempting to create a table with more than INT_MAX cells.");

  }/*THEN*/

  /* allocate and dereference the table. */
  PROTECT(table = allocVector(INTSXP, (int)ncells));
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

/* normalize a conditional probability table. */
SEXP normalize_cpt(SEXP cpt) {

int i = 0, j = 0, nrow = 0, cells = length(cpt);
short int referenced = 0;
double psum = 0;
double *c = NULL;

  /* duplicate the (conditional) probability table if needed... */
  if ((referenced = MAYBE_REFERENCED(cpt)))
    PROTECT(cpt = duplicate(cpt));
  /* ... and dereference it. */
  c = REAL(cpt);

  nrow = INT(getAttrib(cpt, R_DimSymbol));

  for (j = 0; j < (cells/nrow); j++) {

    /* reset column total counter. */
    psum = 0;
    /* compute the new column total. */
    for (i = 0; i < nrow; i++)
      psum += c[CMC(i, j, nrow)];
    /* divide by the new column total. */
    for (i = 0; i < nrow; i++)
      c[CMC(i, j, nrow)] /= psum;

  }/*FOR*/

  if (referenced)
    UNPROTECT(1);

  return cpt;

}/*NORMALIZE_CPT*/

