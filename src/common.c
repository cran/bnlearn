#include "include/rcore.h"
#include "include/matrix.h"

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

int i = 0, *l = NULL;
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

  /* match the elements of the vector against the levels. */
  PROTECT(result = match(levels, vector, 0));

  /* set the levels of the factor. */
  PROTECT(lvls = coerceVector(levels, STRSXP));
  setAttrib(result, R_LevelsSymbol, lvls);

  /* set the class of the return value. */
  setAttrib(result, R_ClassSymbol, mkString("factor"));

  UNPROTECT(3);

  return result;

}/*INT2FAC*/

/* efficient copying of data from the data frame to a matrix for the QR decomposition. */
SEXP qr_matrix(SEXP dataframe, SEXP name) {

SEXP try, result, colnames = getAttrib(dataframe, R_NamesSymbol);
int i = 0, ncols = length(name), nrows = length(VECTOR_ELT(dataframe, 0));
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

int i = 0, j = 0, nrows = 0, cells = length(cpt);
short int duplicated = 0;
double psum = 0;
double *c = NULL;

  /* duplicate the (conditional) probability table if needed... */
  if ((duplicated = NAMED(cpt)) > 0)
    PROTECT(cpt = duplicate(cpt));
  /* ... and dereference it. */
  c = REAL(cpt);

  nrows = INT(getAttrib(cpt, R_DimSymbol));

  for (j = 0; j < (cells/nrows); j++) {

    /* reset column total counter. */
    psum = 0;
    /* compute the new column total. */
    for (i = 0; i < nrows; i++)
      psum += c[CMC(i, j, nrows)];
    /* divide by the new column total. */
    for (i = 0; i < nrows; i++)
      c[CMC(i, j, nrows)] /= psum;

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

  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(obj, R_DimNamesSymbol, dimnames);
  UNPROTECT(1);

}/*SETDIMNAMES*/
