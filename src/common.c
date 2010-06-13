#include "common.h"
#include <R_ext/Arith.h>

/* a rudimental C implementation of which.max(). */
int which_max(double *array, int length) {

int i = 0, imax = 0;
double max = 0;

  for (i = 0; i < length; i++) {

    /* NA and Nan cannot be compared with valid real numbers. */
    if (ISNAN(array[i]))
      return NA_INTEGER;

    if (array[i] > max) {

      imax = i;
      max = array[i];

    }/*THEN*/

  }/*FOR*/

  return imax + 1;

}/*WHICH_MAX*/

/* get the list element named str, or return NULL. */
SEXP getListElement(SEXP list, char *str) {

SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
int i = 0;

  for (i = 0; i < length(list); i++) {

    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {

      elmt = VECTOR_ELT(list, i);
      break;

    }/*THEN*/

  }/*FOR*/

return elmt;

}/*GETLISTELEMENT*/

/* sampling without replacement, internal copy of the
 * SampleNoReplace function in src/main/random.c.
 * Copyright (C) 2003--2008  The R Foundation,
 * licensed under "GPLv2 or later" licence. */
void SampleNoReplace(int k, int n, int *y, int *x) {

int i = 0, j = 0;

  for (i = 0; i < n; i++)
    x[i] = i;
  for (i = 0; i < k; i++) {

    j = n * unif_rand();
    y[i] = x[j] + 1;
    x[j] = x[--n];

  }/*FOR*/

}/*SAMPLENOREPLACE*/

/* return the unique elements from an input vector.*/
SEXP unique(SEXP array) {

int *d = NULL, i = 0, k = 0, dup_counter = 0, n = LENGTH(array);
SEXP dup, result = R_NilValue;

  PROTECT(dup = duplicated(array, FALSE));
  d = LOGICAL(dup);

  for (i = 0; i < n; i++)
    if (d[i] == 0)
      dup_counter++;

  switch(TYPEOF(array)) {

    case INTSXP:

      PROTECT(result = allocVector(INTSXP, dup_counter));
      int *res = INTEGER(result), *a = INTEGER(array);

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          res[k++] = a[i];

      break;

    case STRSXP:

      PROTECT(result = allocVector(STRSXP, dup_counter));

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          SET_STRING_ELT(result, k++, STRING_ELT(array, i));

      break;

    default:

      error("this SEXP type is no handled in unique().");

  }/*SWITCH*/

  UNPROTECT(2);

  return result;

}/*UNIQUE*/

/* determine which elements are dupes. */
SEXP dupe(SEXP array) {

int i = 0, n = LENGTH(array);
int *res = NULL, *tmp = NULL;
SEXP result, temp;

  PROTECT(result = duplicated(array, FALSE));
  PROTECT(temp = duplicated(array, TRUE));
  res = LOGICAL(result);
  tmp = LOGICAL(temp);

  for (i = 0; i < n; i++)
    res[i] = res[i] || tmp[i];

  UNPROTECT(2);

  return result;

}/*DUPE*/

/* transform an integer vector into a factor. */
SEXP int2fac(SEXP vector) {

SEXP result, class, levels, lvls;

  /* find out which levels are needed. */
  PROTECT(levels = unique(vector));

  /* match the elements of the vector against the levels. */
  PROTECT(result = match(levels, vector, 0));

  /* set the levels of the factor. */
  PROTECT(lvls = coerceVector(levels, STRSXP));
  setAttrib(result, R_LevelsSymbol, lvls);

  /* set the class of the return value. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(result, R_ClassSymbol, class);

  UNPROTECT(4);

  return result;

}/*INT2FAC*/

/* in-place conversion an list into a data frame. */
SEXP minimal_data_frame(SEXP obj) {

int n = LENGTH(VECTOR_ELT(obj, 0));
int *row = NULL;
SEXP classname, rownames;

  // generate and set the row names.
  PROTECT(rownames = allocVector(INTSXP, 2));
  row = INTEGER(rownames);

  row[0] = NA_INTEGER;
  row[1] = n;

  setAttrib(obj, R_RowNamesSymbol, rownames);

  // set the class name.      
  PROTECT(classname = allocVector(STRSXP, 1));
  SET_STRING_ELT(classname, 0, mkChar("data.frame"));  
  setAttrib(obj, R_ClassSymbol, classname);

  UNPROTECT(2);

  return obj;

}/*MINIMAL_DATA_FRAME*/

/* efficient column extraction from data frames. */
SEXP dataframe_column(SEXP dataframe, SEXP name, SEXP drop) {

SEXP try, result, colnames = getAttrib(dataframe, R_NamesSymbol);
int *idx = NULL, *d = LOGICAL(drop), nnames = LENGTH(name);


  switch(TYPEOF(name)) {

    case STRSXP:

      /* column names passed as strings; match the corresponding indexes. */
      PROTECT(try = match(colnames, name, 0));
      idx = INTEGER(try);
      break;

    case INTSXP:

      /* these are already indexes, nothing to do. */
      idx = INTEGER(name);
      break;

    default:

      error("this SEXP type is no handled in minimal.data.frame.column().");

  }/*SWITCH*/

  if ((nnames > 1) || (*d == 0)) {

    PROTECT(result = allocVector(VECSXP, nnames));

    for (int i = 0, k = 0; i < nnames; i++)
      SET_VECTOR_ELT(result, k++, VECTOR_ELT(dataframe, idx[i] -1));

    UNPROTECT(1);

  }/*THEN*/
  else {

    result = VECTOR_ELT(dataframe, *idx - 1);

  }/*ELSE*/

  if (TYPEOF(name) != INTSXP)
    UNPROTECT(1);

  return result;

}/*DATAFRAME_COLUMN*/

/* efficient copying of data from the data frame to a matrix for the QR decomposition. */
SEXP qr_matrix(SEXP dataframe, SEXP name) {

SEXP try, result, colnames = getAttrib(dataframe, R_NamesSymbol);
int i = 0, ncols = LENGTH(name), nrows = LENGTH(VECTOR_ELT(dataframe, 0));
int *idx = NULL;
double *res = NULL, *source = NULL;

  PROTECT(try = match(colnames, name, 0));
  idx = INTEGER(try);

  PROTECT(result = allocMatrix(REALSXP, nrows, ncols + 1));
  res = REAL(result);

  /* fill the intercept column. */
  for (i = 0; i < nrows; i++)
    res[i] = 1;

  /* copy the data from the data frame. */
  for (i = 0; i < ncols; i++) {

    /* shift to the beginning of the column. */
    res += nrows;
    /* extract the right column from the data frame. */
    source = REAL(VECTOR_ELT(dataframe, idx[i] -1));
    /* copy the data to the right comun of the matrix. */
    memcpy(res, source, nrows * sizeof(double));

  }/*FOR*/

  UNPROTECT(2);

  return result;

}/*QR_MATRIX*/
