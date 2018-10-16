#include "include/rcore.h"
#include "include/data.frame.h"
#include "include/blas.h"

/* computationally efficient function to compute least squares. */
SEXP fast_lm(SEXP data, SEXP node, SEXP parents, SEXP keep, SEXP missing) {

int i = 0, n = 0, ncol = length(parents);
double **x = NULL, *y = NULL;
SEXP data_x, result, response, coefficients, sd, residuals, fitted, coefnames;

  /* dereference the response variable. */
  PROTECT(response = c_dataframe_column(data, node, TRUE, FALSE));
  y = REAL(response);
  n = length(response);
  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(VECSXP, 4));
  setAttrib(result, R_NamesSymbol,
    mkStringVec(4, "coefficients", "sd", "residuals", "fitted.values"));
  /* allocate the coefficients and the standard error. */
  PROTECT(coefficients = allocVector(REALSXP, ncol + 1));
  PROTECT(coefnames = allocVector(STRSXP, ncol + 1));
  SET_STRING_ELT(coefnames, 0, mkChar("(Intercept)"));
  for (i = 1; i < ncol + 1; i++)
    SET_STRING_ELT(coefnames, i, STRING_ELT(parents, i - 1));
  setAttrib(coefficients, R_NamesSymbol, coefnames);
  PROTECT(sd = allocVector(REALSXP, 1));
  /* extract the relevant columns and dereference them. */
  if (ncol > 0) {

    PROTECT(data_x = c_dataframe_column(data, parents, FALSE, FALSE));
    x = Calloc1D(ncol, sizeof(double *));
    for (i = 0; i < ncol; i++)
      x[i] = REAL(VECTOR_ELT(data_x, i));

  }/*THEN*/

  if (isTRUE(keep)) {

    /* fitted values and residuals are extracted and returned. */
    PROTECT(fitted = allocVector(REALSXP, n));
    PROTECT(residuals = allocVector(REALSXP, n));
    /* estimate all relevant quantities via least squares. */
    c_ols(x, y, n, ncol, REAL(fitted), REAL(residuals),
      REAL(coefficients), REAL(sd), isTRUE(missing));

  }/*THEN*/
  else {

    /* fitted values and residuals are just dummy NAs. */
    fitted = residuals = ScalarReal(NA_REAL);
    /* estimate regression coefficients and standard error via least squares. */
    c_ols(x, y, n, ncol, NULL, NULL, REAL(coefficients), REAL(sd),
      isTRUE(missing));

  }/*ELSE*/

  if (ncol > 0)
    Free1D(x);

  /* save the results and return. */
  SET_VECTOR_ELT(result, 0, coefficients);
  SET_VECTOR_ELT(result, 1, sd);
  SET_VECTOR_ELT(result, 2, residuals);
  SET_VECTOR_ELT(result, 3, fitted);

  UNPROTECT(5 + (ncol > 0) + 2 * isTRUE(keep));

  return result;

}/*FAST_LM*/

/* computationally efficient function to compute conditional least squares. */
SEXP fast_cglm(SEXP data, SEXP node, SEXP parents, SEXP configs, SEXP keep,
    SEXP missing) {

int i = 0, n = 0, *z = NULL, ncol = length(parents), nz = 0;
double **x = NULL, *y = NULL;
SEXP response, coefficients, sd, residuals, fitted, coefnames;
SEXP data_x, result, confnames, dummy_configs;

  /* dereference the response variable and the configurations. */
  PROTECT(response = c_dataframe_column(data, node, TRUE, FALSE));
  y = REAL(response);
  n = length(response);
  z = INTEGER(configs);
  /* prepare the configuration index names. */
  confnames = getAttrib(configs, R_LevelsSymbol);
  nz = length(confnames);
  /* allocate the coefficients and the standard errors. */
  PROTECT(coefficients = allocMatrix(REALSXP, ncol + 1, nz));
  PROTECT(coefnames = allocVector(STRSXP, ncol + 1));
  SET_STRING_ELT(coefnames, 0, mkChar("(Intercept)"));
  for (i = 1; i < ncol + 1; i++)
    SET_STRING_ELT(coefnames, i, STRING_ELT(parents, i - 1));
  setDimNames(coefficients, coefnames, confnames);
  PROTECT(sd = allocVector(REALSXP, nz));
  setAttrib(sd, R_NamesSymbol, confnames);
  /* extract the relevant columns and dereference them. */
  if (ncol > 0) {

    PROTECT(data_x = c_dataframe_column(data, parents, FALSE, FALSE));
    x = Calloc1D(ncol, sizeof(double *));
    for (i = 0; i < ncol; i++)
      x[i] = REAL(VECTOR_ELT(data_x, i));

  }/*THEN*/

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(VECSXP, 5));
  setAttrib(result, R_NamesSymbol,
    mkStringVec(5, "coefficients", "sd", "configs", "residuals", "fitted.values"));
  if (isTRUE(keep)) {

    /* configurations are copied from the data. */
    SET_VECTOR_ELT(result, 2, configs);
    /* fitted values and residuals are extracted and returned. */
    PROTECT(fitted = allocVector(REALSXP, n));
    PROTECT(residuals = allocVector(REALSXP, n));
    /* estimate all relevant quantities via least squares. */
    c_cls(x, y, z, n, ncol, nz, REAL(fitted), REAL(residuals),
      REAL(coefficients), REAL(sd), isTRUE(missing));

  }/*THEN*/
  else {

    /* fitted values, residuals and configurations are just dummy NAs. */
    PROTECT(fitted = ScalarReal(NA_REAL));
    PROTECT(residuals = ScalarReal(NA_REAL));
    PROTECT(dummy_configs = allocVector(INTSXP, 1));
    INT(dummy_configs) = NA_INTEGER;
    setAttrib(dummy_configs, R_ClassSymbol, mkString("factor"));
    setAttrib(dummy_configs, R_LevelsSymbol, confnames);
    SET_VECTOR_ELT(result, 2, dummy_configs);
    /* estimate all relevant quantities via least squares. */
    c_cls(x, y, z, n, ncol, nz, NULL, NULL, REAL(coefficients), REAL(sd),
      isTRUE(missing));

  }/*ELSE*/

  if (ncol > 0)
    Free1D(x);

  /* save the results and return. */
  SET_VECTOR_ELT(result, 0, coefficients);
  SET_VECTOR_ELT(result, 1, sd);
  SET_VECTOR_ELT(result, 3, residuals);
  SET_VECTOR_ELT(result, 4, fitted);

  UNPROTECT(5 + (ncol > 0) + (isTRUE(keep) ? 2 : 3));

  return result;

}/*FAST_CGLM*/

