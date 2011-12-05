#include "common.h"
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>

/* a rudimental C implementation of which.max(). */
int which_max(double *array, int length) {

int i = 0, imax = -1;
double max = R_NegInf;

  for (i = 0; i < length; i++) {

    /* NA and NaN cannot be compared with valid real numbers. */
    if (ISNAN(array[i]))
      continue;

    if (array[i] > max) {

      imax = i;
      max = array[i];

    }/*THEN*/

  }/*FOR*/

    /* if all elements are NA/NaN return NA. */
    if (imax < 0)
      return NA_INTEGER;

  return imax + 1;

}/*WHICH_MAX*/

/* return all maxima in the array, modulo numeric tolerance. */
void all_max(double *array, int length, int *maxima, int *nmax, int *indexes) {

  int i = 0;
  double tol = MACHINE_TOL;

  /* sort the elements of the array. */
  rsort_with_index(array, indexes, length);

  /* count the number of maxima (considering numeric tolerance). */
  for (i = length - 1; i >= 0; i--)
    if (array[i] < array[length - 1] - tol)
      break;

  /* set the counter for the number of maxima. */
  *nmax = length - i - 1;

  /* save the indexes of the maxima. */
  memcpy(maxima, indexes + length - *nmax, *nmax * sizeof(int));

}/*ALL_MAX*/

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

      error("this SEXP type is not handled in unique().");

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
SEXP int2fac(SEXP vector, int *nlevels) {

int i = 0, *l = NULL;
SEXP result, class, levels, lvls;

  if (!nlevels) {

    PROTECT(levels = unique(vector));

  }/*THEN*/
  else {

    PROTECT(levels = allocVector(INTSXP, *nlevels));
    l = INTEGER(levels);

    for (i = 0; i < *nlevels; i++)
      l[i] = i;

  }/*ELSE*/

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
int *idx = NULL, *d = LOGICAL(drop);
int nnames = LENGTH(name), name_type = TYPEOF(name);

  if (dataframe == R_NilValue)
    return R_NilValue;

  switch(name_type) {

    case STRSXP:

      /* column names passed as strings; match the corresponding indexes. */
      PROTECT(try = match(colnames, name, 0));
      idx = INTEGER(try);
      break;

    case REALSXP:

      /* these are almost good enough, coerce them to integers. */
      PROTECT(try = coerceVector(name, INTSXP));
      idx = INTEGER(try);
      break;

    case INTSXP:

      /* these are already indexes, nothing to do. */
      idx = INTEGER(name);
      break;

    default:

      error("this SEXP type is not handled in minimal.data.frame.column().");

  }/*SWITCH*/

  if ((nnames > 1) || (*d == 0)) {

    PROTECT(result = allocVector(VECSXP, nnames));

    for (int i = 0, k = 0; i < nnames; i++)
      SET_VECTOR_ELT(result, k++, VECTOR_ELT(dataframe, idx[i] -1));

    UNPROTECT(1);

  }/*THEN*/
  else {

    if (*idx != 0)
      result = VECTOR_ELT(dataframe, *idx - 1);
    else
      result = R_NilValue;

  }/*ELSE*/

  if (name_type != INTSXP)
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

/* set the column names on arc sets. */
SEXP finalize_arcs(SEXP arcs) {

SEXP dimnames, colnames;

  /* allocate column names. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* set the column names. */
  setAttrib(arcs, R_DimNamesSymbol, dimnames);

  UNPROTECT(2);

  return arcs;

}/*FINALIZE_ARCS*/

/* inverse function of the UPTRI3() macro. */
void INV_UPTRI3(int x, int n, int *res) {

int c = 0, r = 0, cn = n - 1;

  for (r = 0; r < n; r++) {

    if (x < cn) {
      c = n - (cn - x);
      break;
    }
    else
      cn += n - (r + 2);

  }/*FOR*/

  res[0] = r;
  res[1] = c;

}/*INV_UPTRI3*/

