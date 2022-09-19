#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../math/linear.algebra.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../minimal/common.h"

SEXP gaussian_ols_parameters(SEXP data, SEXP node, SEXP parents, SEXP keep,
    SEXP replace_unidentifiable, SEXP missing) {

int i = 0, n = 0, ncol = length(parents);
double **x = NULL, *y = NULL, *coefs = NULL;
double sd = 0;
SEXP data_x, result, response, coefficients, residuals, fitted, coefnames;

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
  coefs = REAL(coefficients);
  PROTECT(coefnames = allocVector(STRSXP, ncol + 1));
  SET_STRING_ELT(coefnames, 0, mkChar("(Intercept)"));
  for (i = 1; i < ncol + 1; i++)
    SET_STRING_ELT(coefnames, i, STRING_ELT(parents, i - 1));
  setAttrib(coefficients, R_NamesSymbol, coefnames);
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
      coefs, &sd, isTRUE(missing));

  }/*THEN*/
  else {

    /* fitted values and residuals are just dummy NAs. */
    fitted = residuals = ScalarReal(NA_REAL);
    /* estimate regression coefficients and standard error via least squares. */
    c_ols(x, y, n, ncol, NULL, NULL, coefs, &sd, isTRUE(missing));

  }/*ELSE*/

  /* replace unidentifiable regression coefficients and the standard error with
   * zeroes to prevent NAs from propagating. */
  if (isTRUE(replace_unidentifiable)) {

    for (i = 0; i < ncol + 1; i++)
      if (ISNAN(coefs[i]))
        coefs[i] = 0;

    if (ISNAN(sd))
      sd = 0;

  }/*THEN*/

  if (ncol > 0)
    Free1D(x);

  /* save the results and return. */
  SET_VECTOR_ELT(result, 0, coefficients);
  SET_VECTOR_ELT(result, 1, mkReal(sd));
  SET_VECTOR_ELT(result, 2, residuals);
  SET_VECTOR_ELT(result, 3, fitted);

  UNPROTECT(4 + (ncol > 0) + 2 * isTRUE(keep));

  return result;

}/*GAUSSIAN_OLS_PARAMETERS*/
