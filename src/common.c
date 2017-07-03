#include "include/rcore.h"
#include "include/matrix.h"
#include "include/sets.h"

/* class check to match that in the R code. */
int c_is(SEXP obj, const char *str) {

int i = 0;
SEXP class;

  class = getAttrib(obj, R_ClassSymbol);

  for (i = 0; i < length(class); i++)
    if (!strcmp(CHAR(STRING_ELT(class, i)), str))
      return TRUE;

  return FALSE;

}/*C_IS*/

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

int *d = NULL, i = 0, k = 0, dup_counter = 0, n = length(array);
int *res = NULL, *a = NULL;
SEXP dup, result = R_NilValue;

  PROTECT(dup = duplicated(array, FALSE));
  d = LOGICAL(dup);

  switch(TYPEOF(array)) {

    case INTSXP:

      a = INTEGER(array);

      for (i = 0; i < n; i++)
        if ((d[i] == 0) && (a[i] != NA_INTEGER))
          dup_counter++;

      PROTECT(result = allocVector(INTSXP, dup_counter));
      res = INTEGER(result);

      for (i = 0; i < n; i++)
        if ((d[i] == 0) && (a[i] != NA_INTEGER))
          res[k++] = a[i];

      break;

    case STRSXP:

      for (i = 0; i < n; i++)
        if (d[i] == 0)
          dup_counter++;

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

int i = 0, n = length(array);
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

int i = 0, *l = NULL, *r = NULL, *v = INTEGER(vector);
SEXP result, levels, lvls;

  if (!nlevels) {

    PROTECT(levels = unique(vector));

  }/*THEN*/
  else {

    PROTECT(levels = allocVector(INTSXP, *nlevels));
    l = INTEGER(levels);

    for (i = 0; i < *nlevels; i++)
      l[i] = i;

  }/*ELSE*/

  /* match the elements of the vector against the levels, preserving NAs and
   * setting to NA those elements that cannot be matched to a level. */
  PROTECT(result = match(levels, vector, 0));
  r = INTEGER(result);

  for (i = 0; i < length(result); i++)
    if ((r[i] == 0) || (v[i] == NA_INTEGER))
      r[i] = NA_INTEGER;

  /* set the levels of the factor. */
  PROTECT(lvls = coerceVector(levels, STRSXP));
  setAttrib(result, R_LevelsSymbol, lvls);

  /* set the class of the return value. */
  setAttrib(result, R_ClassSymbol, mkString("factor"));

  UNPROTECT(3);

  return result;

}/*INT2FAC*/

/* efficient copying of data from the data frame to a matrix for the QR decomposition. */
void c_qr_matrix(double *qr, double **x, int nrow, int ncol) {

int i = 0;

  /* fill the intercept column. */
  for (i = 0; i < nrow; i++)
    qr[i] = 1;

  /* copy the data to the right comun of the matrix. */
  for (i = 0; i < ncol; i++)
    memcpy(qr + (i + 1) * nrow, x[i], nrow * sizeof(double));

}/*C_QR_MATRIX*/

/* inverse function of the UPTRI3() macro. */
void INV_UPTRI3(int x, int n, int *res) {

int c = 0, r = 0, cn = n - 1;

  for (r = 0; r < n; r++) {

    if (x < cn) {

      c = n - (cn - x);
      break;

    }/*THEN*/
    else {

      cn += n - (r + 2);

    }/*ELSE*/

  }/*FOR*/

  res[0] = r;
  res[1] = c;

}/*INV_UPTRI3*/

/* normalize a conditional probability table. */
SEXP normalize_cpt(SEXP cpt) {

int i = 0, j = 0, nrow = 0, cells = length(cpt);
short int duplicated = 0;
double psum = 0;
double *c = NULL;

  /* duplicate the (conditional) probability table if needed... */
  if ((duplicated = NAMED(cpt)) > 0)
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

  if (duplicated > 0)
    UNPROTECT(1);

  return cpt;

}/*NORMALIZE_CPT*/

/* data-independent pointer extraction without USE_RINTERNALS. */
void *DATAPTR(SEXP x) {

  switch(TYPEOF(x)) {

    case REALSXP:

      return (void *)REAL(x);

    case INTSXP:

      return (void *)INTEGER(x);

  }/*SWITCH*/

  return NULL;

}/*DATAPTR*/

/* variadic version of mkString(). */
SEXP mkRealVec(int n, ...) {

va_list strings;
int i = 0;
double *v = NULL;
SEXP vec;

  PROTECT(vec = allocVector(REALSXP, n));
  v = REAL(vec);
  va_start(strings, n);
  for (i = 0; i < n; i++)
    v[i] = va_arg(strings, double);
  va_end(strings);
  UNPROTECT(1);

  return vec;

}/*MKSTRINGVEC*/

/* set row and column names. */
void setDimNames(SEXP obj, SEXP rownames, SEXP colnames) {

SEXP dimnames;

  PROTECT(rownames);
  PROTECT(colnames);
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(obj, R_DimNamesSymbol, dimnames);
  UNPROTECT(3);

}/*SETDIMNAMES*/

/* minimal implementation of table(). */
SEXP minimal_table(SEXP dataframe) {

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

  for (i = 0; i < nrow; i++)
    tt[cfg[i]]++;

  /* set the attributess for class and dimensions. */
  setAttrib(table, R_ClassSymbol, mkString("table"));
  setAttrib(table, R_DimSymbol, dims);
  setAttrib(table, R_DimNamesSymbol, dimnames);

  UNPROTECT(3);

  Free1D(columns);
  Free1D(cfg);

  return table;

}/*MINIMAL_TABLE*/
