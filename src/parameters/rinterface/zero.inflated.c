#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../include/globals.h"
#include "../../math/hyperpoisson.h"
#include "../../math/xnegbin.h"
#include "../../minimal/common.h"
#include "../../minimal/data.frame.h"
#include "../../minimal/strings.h"
#include "../../scores/scores.h"

#define coefs_to_params(family, zprob, zcoefs, counts, ccoefs, data) \
  do { \
  if (family == GLM_HYPERPOISSON) \
    zihp_coefs_to_params(zprob, zcoefs, counts, ccoefs, data); \
  else \
    zinb_coefs_to_params(zprob, zcoefs, counts, ccoefs, data); \
  } while (0)

/* package the fitted parameters of a zero-inflated count node to pass to R. */
static SEXP package_zi_result(glm_family_e family, double *estimates, SEXP node,
    SEXP parents, SEXP data, SEXP keep, SEXP replace_unidentifiable,
    int nparents, int nobs) {

double *zin = NULL, *cn2 = NULL, scalar = 0;
double *fit = NULL, *resid = NULL, *zinf_prob = NULL, *canon = NULL, *y = NULL;
const char *comp2name = NULL, *scalname = NULL;
bool saturated = FALSE;
SEXP result, zeroinfl, comp2, coefnames, fitted, residuals, par_data;
tabular par_df = { 0 };


  /* set the parameter set names in the return value. */
  if (family == GLM_HYPERPOISSON) {

    comp2name = "intensity";
    scalname = "dispersion";

  }/*THEN*/
  else if (family == GLM_NEGBIN) {

    comp2name = "prsucc";
    scalname = "failures";

  }/*ELSE*/

  /* allocate the return value. */
  PROTECT(result = allocVector(VECSXP, 5));
  setAttrib(result, R_NamesSymbol,
    mkStringVec(5, "inflation", comp2name, scalname, "residuals",
                   "fitted.values"));
  PROTECT(zeroinfl = allocVector(REALSXP, nparents + 1));
  zin = REAL(zeroinfl);
  PROTECT(comp2 = allocVector(REALSXP, nparents + 1));
  cn2 = REAL(comp2);

  /* allocate the coefficient names and set them. */
  PROTECT(coefnames = allocVector(STRSXP, nparents + 1));
  SET_STRING_ELT(coefnames, 0, mkChar("(Intercept)"));
  for (int i = 1; i < nparents + 1; i++)
    SET_STRING_ELT(coefnames, i, STRING_ELT(parents, i - 1));
  setAttrib(zeroinfl, R_NamesSymbol, coefnames);
  setAttrib(comp2, R_NamesSymbol, coefnames);

  /* copy the parameters over. */
  memcpy(zin, estimates, (nparents + 1) * sizeof(double));
  memcpy(cn2, estimates + nparents + 1, (nparents + 1) * sizeof(double));
  scalar = exp(estimates[2 * nparents + 2]);

  if (isTRUE(keep)) {

    /* extract and set the fitted values and residuals. */
    PROTECT(fitted = allocVector(REALSXP, nobs));
    PROTECT(residuals = allocVector(REALSXP, nobs));
    fit = REAL(fitted);
    resid = REAL(residuals);

    PROTECT(par_data = c_dataframe_column(data, parents, FALSE, TRUE));
    par_df = tabular_from_SEXP(par_data, 0, 0);

    zinf_prob = Calloc1D(nobs, sizeof(double));
    canon = Calloc1D(nobs, sizeof(double));

    /* translate the regression coefficients into canonical parameters. */
    coefs_to_params(family, zinf_prob, zin, canon, cn2, par_df);

    if (family == GLM_HYPERPOISSON) {

      for (int i = 0; i < nobs; i++)
        fit[i] = (1 - zinf_prob[i]) * hypois_mean(canon[i], scalar);

    }/*THEN*/
    else {

      for (int i = 0; i < nobs; i++)
        fit[i] = (1 - zinf_prob[i]) * negbin_mean(canon[i], scalar);

    }/*ELSE*/

    y = REAL(c_dataframe_column(data, node, TRUE, FALSE));

    for (int i = 0; i < nobs; i++)
      resid[i] = y[i] - fit[i];

    FreeTAB(par_df);
    Free1D(zinf_prob);
    Free1D(canon);

    UNPROTECT(1);

  }/*THEN*/
  else {

    /* fitted values and residuals are just dummy NAs. */
    fitted = residuals = ScalarReal(NA_REAL);

  }/*ELSE*/

  /* detect saturated fits, but leave the other boundaries (vanishing intensity
   * or success probability, the inflation probability close to 0/1) alone. */
  if (isTRUE(replace_unidentifiable)) {

    double *zp = Calloc1D(nobs, sizeof(double));
    double *cp = Calloc1D(nobs, sizeof(double));
    SEXP sat_data;
    tabular sat_df = { 0 };

    PROTECT(sat_data = c_dataframe_column(data, parents, FALSE, TRUE));
    sat_df = tabular_from_SEXP(sat_data, 0, 0);

    coefs_to_params(family, zp, zin, cp, cn2, sat_df);

    if (family == GLM_HYPERPOISSON) {

      for (int i = 0; i < nobs; i++)
        if (!R_FINITE(cp[i]) && !ISNAN(cp[i])) {
          saturated = TRUE;
          break;
        }/*THEN*/

    }/*THEN*/
    else if (family == GLM_NEGBIN) {

      for (int i = 0; i < nobs; i++)
        if (cp[i] >= 1 - 10 * MACHINE_TOL) {
          saturated = TRUE;
          break;
        }/*THEN*/

    }/*THEN*/

    FreeTAB(sat_df);
    Free1D(zp);
    Free1D(cp);
    UNPROTECT(1);

  }/*THEN*/

  /* replace unidentifiable regression coefficients and the standard error with
   * zeroes to prevent NAs and saturated estimates from propagating. */
  if (isTRUE(replace_unidentifiable)) {

    for (int i = 0; i < nparents + 1; i++)
      if (ISNAN(zin[i]))
        zin[i] = 0;
    for (int i = 0; i < nparents + 1; i++)
      if (ISNAN(cn2[i]) || saturated)
        cn2[i] = 0;
    if (!R_FINITE(scalar) || saturated)
      scalar = 1;

  }/*THEN*/

  /* save the results and return. */
  SET_VECTOR_ELT(result, 0, zeroinfl);
  SET_VECTOR_ELT(result, 1, comp2);
  SET_VECTOR_ELT(result, 2, mkReal(scalar));
  SET_VECTOR_ELT(result, 3, residuals);
  SET_VECTOR_ELT(result, 4, fitted);

  UNPROTECT(4 + 2 * isTRUE(keep));

  return result;

}/*PACKAGE_ZI_RESULT*/

SEXP zero_inflated_hyperpoisson_parameters(SEXP node, SEXP parents, SEXP dag,
    SEXP data, SEXP keep, SEXP replace_unidentifiable, SEXP missing,
    SEXP em_max_iter, SEXP em_tol, SEXP one_step) {

int nparents = length(parents), nobs = length(VECTOR_ELT(data, 0));
double *estimates = NULL;
SEXP result;

  estimates = Calloc1D(2 * nparents + 3, sizeof(double));

  /* fit the parameters by EM with the GLM components fitted by IRLS. */
  em_irls_node(node, dag, data, estimates, NULL, isTRUE(missing), 0, FALSE,
    asInteger(em_max_iter), asReal(em_tol), isTRUE(one_step), GLM_HYPERPOISSON);

  result = package_zi_result(GLM_HYPERPOISSON, estimates, node, parents, data,
             keep, replace_unidentifiable, nparents, nobs);

  Free1D(estimates);

  return result;

}/*ZERO_INFLATED_HYPERPOISSON_PARAMETERS*/

SEXP zero_inflated_negative_binomial_parameters(SEXP node, SEXP parents, SEXP dag,
    SEXP data, SEXP keep, SEXP replace_unidentifiable, SEXP missing,
    SEXP em_max_iter, SEXP em_tol, SEXP one_step) {

int nparents = length(parents), nobs = length(VECTOR_ELT(data, 0));
double *estimates = NULL;
SEXP result;

  estimates = Calloc1D(2 * nparents + 3, sizeof(double));

  /* fit the parameters by EM with the GLM components fitted by IRLS. */
  em_irls_node(node, dag, data, estimates, NULL, isTRUE(missing), 0, FALSE,
    asInteger(em_max_iter), asReal(em_tol), isTRUE(one_step), GLM_NEGBIN);

  result = package_zi_result(GLM_NEGBIN, estimates, node, parents, data,
             keep, replace_unidentifiable, nparents, nobs);

  Free1D(estimates);

  return result;

}/*ZERO_INFLATED_NEGATIVE_BINOMIAL_PARAMETERS*/
