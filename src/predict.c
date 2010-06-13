#include "common.h"

/* predict the value of a gaussian node without parents. */
SEXP gpred(SEXP fitted, SEXP data) {

int i = 0, *ndata = INTEGER(data);
double *mean = NULL, *res = NULL;
SEXP result;

  /* get the (only) coefficient of the linear regression. */
  mean = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, *ndata));
  res = REAL(result);

  /* copy the mean in the return value. */
  for (i = 0; i < *ndata; i++)
    res[i] = *mean;

  UNPROTECT(1);

  return result;

}/*GPRED*/

/* predict the value of a gaussian node with one or more parents. */
SEXP cgpred(SEXP fitted, SEXP data)  {

int i = 0, j = 0, ndata = LENGTH(VECTOR_ELT(data, 0)), ncols = LENGTH(data);
double *res = NULL, *coefs = NULL;
double **columns = NULL;
SEXP result;

  /* get the coefficient of the linear regression. */
  coefs = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, ndata));
  res = REAL(result);

  /* dereference the columns of the data frame. */
  columns = Calloc(ncols, double *);
  for (i = 0; i < ncols; i++)
    columns[i] = REAL(VECTOR_ELT(data, i));

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    res[i] = coefs[0];

    for (j = 0; j < ncols; j++)
      res[i] += columns[j][i] * coefs[j + 1];

  }/*FOR*/

  Free(columns);

  UNPROTECT(1);

  return result;

}/*CGPRED*/

/* predict the value of a discrete node without parents. */
SEXP dpred(SEXP fitted, SEXP data) {

int i = 0, imax = 0, ndata = LENGTH(data);
int *res = NULL;
double *prob = NULL; 
SEXP ptab, result;

  /* get the probabilities of the multinomial distribution. */
  ptab = getListElement(fitted, "prob");
  prob = REAL(ptab);

  /* find out the mode. */
  imax = which_max(prob, LENGTH(ptab)) + 1;

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  /* copy the index of the mode in the return value. */
  for (i = 0; i < ndata; i++)
    res[i] = imax;

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, getAttrib(data, R_LevelsSymbol));
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*DPRED*/

/* predict the value of a discrete node with one or more parents. */
SEXP cdpred(SEXP fitted, SEXP data, SEXP parents) {

int i = 0, ndata = LENGTH(data), nrows = 0, ncols = 0;
int *configs = INTEGER(parents), *mode = NULL, *res = NULL;
double *prob = NULL;
SEXP temp, result;

  /* get the probabilities of the multinomial distribution. */
  temp = getListElement(fitted, "prob");
  nrows = INT(getAttrib(temp, R_DimSymbol));
  ncols = LENGTH(temp) / nrows;
  prob = REAL(temp);

  /* get the mode for each configuration. */
  mode = alloc1dcont(ncols);
  for (i = 0; i < ncols; i++)
    mode[i] = which_max(prob + nrows * i, nrows);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  /* copy the index of the mode in the return value. */
  for (i = 0; i < ndata; i++)
    res[i] = mode[configs[i] - 1];

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, getAttrib(data, R_LevelsSymbol));
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*CDPRED*/

