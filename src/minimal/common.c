#include "../include/rcore.h"
#include "../minimal/unique.h"
#include "../math/linear.algebra.h"

/* class check to match that in the R code. */
bool c_is(SEXP obj, const char *str) {

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

/* efficient copying of data to a matrix suitable for a QR decomposition. */
void c_qr_matrix(double *qr, double **x, int nrow, int ncol, int *complete,
    int ncomplete) {

int i = 0, j = 0, k = 0;

  if (complete) {

    /* fill the intercept column. */
    for (i = 0; i < ncomplete; i++)
      qr[i] = 1;

    /* copy the data to the right comun of the matrix. */
    for (j = 0; j < ncol; j++)
      for (i = 0, k = 0; i < nrow; i++)
        if (complete[i])
          qr[k++ + (j + 1) * ncomplete] = x[j][i];

  }/*THEN*/
  else {

    /* fill the intercept column. */
    for (i = 0; i < nrow; i++)
      qr[i] = 1;

    /* copy the data to the right comun of the matrix. */
    for (i = 0; i < ncol; i++)
      memcpy(qr + (i + 1) * nrow, x[i], nrow * sizeof(double));

  }/*ELSE*/

}/*C_QR_MATRIX*/

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

/* encapsulate a floating point number in a SEXP. */
SEXP mkReal(double x) {

SEXP scalar;

  PROTECT(scalar = allocVector(REALSXP, 1));
  NUM(scalar) = x;
  UNPROTECT(1);

  return scalar;

}/*MKREAL*/

/* variadic version of mkReal(). */
SEXP mkRealVec(int n, ...) {

va_list reals;
int i = 0;
double *v = NULL;
SEXP vec;

  PROTECT(vec = allocVector(REALSXP, n));
  v = REAL(vec);
  va_start(reals, n);
  for (i = 0; i < n; i++)
    v[i] = va_arg(reals, double);
  va_end(reals);
  UNPROTECT(1);

  return vec;

}/*MKREALVEC*/

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

/* subset a vector using names or integer values. */
SEXP subset_by_name(SEXP vec, int n, ...) {

int i = 0, j = 0, k = 0, idx = 0, subvec_len = 0, cur_type = 0;
va_list strings;
SEXP subvec, subvec_names, cur, try, vec_names;

  /* compute the length of the return value. */
  va_start(strings, n);
  for (i = 0; i < n; i++)
    subvec_len += length(va_arg(strings, SEXP));
  va_end(strings);

  /* allocate the return value. */
  PROTECT(subvec = allocVector(TYPEOF(vec), subvec_len));
  PROTECT(subvec_names = allocVector(STRSXP, subvec_len));
  setAttrib(subvec, R_NamesSymbol, subvec_names);

  PROTECT(vec_names = getAttrib(vec, R_NamesSymbol));
  va_start(strings, n);
  for (i = 0, k = 0; i < n; i++) {

    /* find out which elements to extract... */
    cur = va_arg(strings, SEXP);

    if (isNull(cur))
      continue;

    cur_type = TYPEOF(cur);
    if (cur_type == STRSXP)
      PROTECT(try = match(vec_names, cur, 0));
    else if (cur_type == INTSXP)
      try = cur;
    else
      error("unknown subset object type (class: %s).",
        CHAR(STRING_ELT(getAttrib(cur, R_ClassSymbol), 0)));

    /* ... and save them in the return value. */
    switch(TYPEOF(vec)) {

      case LGLSXP:
        for (j = 0; j < length(cur); j++) {

          idx = INTEGER(try)[j] - 1;
          LOGICAL(subvec)[k] = LOGICAL(vec)[idx];
          if (cur_type == STRSXP)
            SET_STRING_ELT(subvec_names, k++, STRING_ELT(cur, j));
          else if (cur_type == INTSXP)
            SET_STRING_ELT(subvec_names, k++, STRING_ELT(vec_names, idx));

        }/*FOR*/
        break;

      case REALSXP:
        for (j = 0; j < length(cur); j++) {
          idx = INTEGER(try)[j] - 1;
          REAL(subvec)[k] = REAL(vec)[idx];
          if (cur_type == STRSXP)
            SET_STRING_ELT(subvec_names, k++, STRING_ELT(cur, j));
          else if (cur_type == INTSXP)
            SET_STRING_ELT(subvec_names, k++, STRING_ELT(vec_names, idx));

        }/*FOR*/
        break;

    }/*SWITCH*/

    if (cur_type == STRSXP)
      UNPROTECT(1);

  }/*FOR*/
  va_end(strings);

  UNPROTECT(3);

  return subvec;

}/*SUBSET_BY_NAME*/

/* check all elements of a SEXP are equal to a value. */
bool all_equal(SEXP vec, SEXP val) {

int i = 0, *vl = NULL, al = FALSE;

  switch(TYPEOF(vec)) {

    case LGLSXP:
      vl = LOGICAL(vec);
      al = isTRUE(val);
      for (i = 0; i < length(vec); i++)
        if (vl[i] != al)
          return FALSE;
      break;

    default:
      error("unknown object type (class: %s).",
        CHAR(STRING_ELT(getAttrib(vec, R_ClassSymbol), 0)));

  }/*SWITCH*/

  return TRUE;

}/*ALL_EQUAL*/
