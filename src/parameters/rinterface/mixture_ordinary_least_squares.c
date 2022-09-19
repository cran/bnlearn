#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../minimal/common.h"
#include "../../math/linear.algebra.h"

SEXP mixture_gaussian_ols_parameters(SEXP data, SEXP node, SEXP parents,
    SEXP configs, SEXP keep, SEXP replace_unidentifiable, SEXP missing) {

int i = 0, n = 0, *z = NULL, ncol = length(parents), nz = 0;
double **x = NULL, *y = NULL, *coefs = NULL, *sds = NULL;
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
  coefs = REAL(coefficients);
  PROTECT(coefnames = allocVector(STRSXP, ncol + 1));
  SET_STRING_ELT(coefnames, 0, mkChar("(Intercept)"));
  for (i = 1; i < ncol + 1; i++)
    SET_STRING_ELT(coefnames, i, STRING_ELT(parents, i - 1));
  setDimNames(coefficients, coefnames, confnames);
  PROTECT(sd = allocVector(REALSXP, nz));
  setAttrib(sd, R_NamesSymbol, confnames);
  sds = REAL(sd);
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
    c_cls(x, y, z, n, ncol, nz, REAL(fitted), REAL(residuals), coefs, sds,
      isTRUE(missing));

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
    c_cls(x, y, z, n, ncol, nz, NULL, NULL, coefs, sds, isTRUE(missing));

  }/*ELSE*/

  /* replace unidentifiable regression coefficients and the standard error with
   * zeroes to prevent NAs from propagating. */
  if (isTRUE(replace_unidentifiable)) {

    for (i = 0; i < (ncol + 1) * nz; i++)
      if (ISNAN(coefs[i]))
        coefs[i] = 0;

    for (i = 0; i < nz; i++)
      if (ISNAN(sds[i]))
        sds[i] = 0;

  }/*THEN*/

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

