#include "common.h"

/* macro for the [i,j] element of the data frame. */
#define DATAFRAME(i, j) INTEGER(VECTOR_ELT(parents, j))[i]
/* macro for the number of levels of the [,j] column. */
#define NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(parents, j), R_LevelsSymbol))

SEXP cfg(SEXP parents);

/* wrapper around the cfg() function, which returns a factor instead of
 * an integer vector. */
SEXP cfg2(SEXP parents) {

SEXP temp, result;

  PROTECT(temp = cfg(parents));
  result = int2fac(temp);

  UNPROTECT(1);

  return result;

}/*CFG2*/

/* identify different configurations of factors and assign them unique integer
 * codes, computed using the general formula for the column-mayor indexing. */
SEXP cfg(SEXP parents) {

int i = 0, j = 0, cfgmap = 0;
int ncols = LENGTH(parents), nrows = LENGTH(VECTOR_ELT(parents, 0));
int *cumlevels = NULL, *res = NULL;
SEXP result;

  /* create the cumulative products of the number of levels. */
  cumlevels = alloc1dcont(ncols);

  /* set the first one to 1 ... */
  cumlevels[0] = 1;

  /* ... then compute the following ones. */
  for (j = 1; j < ncols; j++)
    cumlevels[j] = cumlevels[j - 1] * NLEVELS(j - 1);

  /* allocate an array of size nrow for the configuration. */
  PROTECT(result = allocVector(INTSXP, nrows));
  res = INTEGER(result);

  for (i = 0; i < nrows; i++) {

    /* reset the configuration mapping of the new row. */
    cfgmap = 0;

    for (j = 0; j < ncols; j++) {

      cfgmap += (DATAFRAME(i, j) - 1) * cumlevels[j];

    }/*FOR*/

  /* save the configuration in the array. */
  res[i] = cfgmap;

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*CFG*/

/* collapse a multidimensional table into a bidimensional table with one column
 * for each configuration of the 2nd, 3rd, etc. variables. */
SEXP collapse_table(SEXP table) {

int i = 0, ncells = 0, nrows = 0, ncols = 0;
int *lab = NULL;
SEXP result, dimnames, rownames, colnames, labels, dimlabels;

  /* compute the dimensions of the table. */
  ncells = LENGTH(table);
  nrows = INT(getAttrib(table, R_DimSymbol));
  ncols = ncells / nrows;

  /* allocate the return value. */
  PROTECT(result = allocMatrix(REALSXP, nrows, ncols));

  /* copy the conditional probabilities from the old to the new table. */
  memcpy(REAL(result), REAL(table), ncells * sizeof(double));

  /* allocate and initialize the row and column names. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  /* row names are the same as in the old table. */
  rownames = VECTOR_ELT(getAttrib(table, R_DimNamesSymbol), 0);
  /* column names are integer numbers converted to strings. */
  PROTECT(labels = allocVector(INTSXP, ncols));
  lab = INTEGER(labels);

  for(i = 0; i < ncols; i++)
    lab[i] = i;

  /* set row and column names. */
  PROTECT(colnames = coerceVector(labels, STRSXP));
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_VECTOR_ELT(dimnames, 1, colnames);
  PROTECT(dimlabels = allocVector(STRSXP, 2));
  SET_STRING_ELT(dimlabels, 0, mkChar("data"));
  SET_STRING_ELT(dimlabels, 1, mkChar("cfg"));
  setAttrib(dimnames, R_NamesSymbol, dimlabels);  
  setAttrib(result, R_DimNamesSymbol, dimnames);  

  UNPROTECT(5);

  return result;

}/*COLLAPSE_TABLE*/
