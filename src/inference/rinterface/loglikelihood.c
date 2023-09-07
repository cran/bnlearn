#include "../../include/rcore.h"
#include "../../include/globals.h"
#include "../../fitted/fitted.h"
#include "../../core/data.table.h"
#include "../../core/allocations.h"
#include "../../minimal/common.h"
#include "../loglikelihood/loglikelihood.h"
#include "../loss.h"

/* multi-purpose log-likelihood function for data and a fitted network. */
SEXP loglikelihood_function(SEXP fitted, SEXP data, SEXP by_sample,
    SEXP keep_nodes, SEXP propagate_missing, SEXP debug) {

int ndata = length(VECTOR_ELT(data, 0));
double *loglik = NULL;
fitted_bn bn = fitted_network_from_SEXP(fitted);
bool by = isTRUE(by_sample), propagate = isTRUE(propagate_missing);
bool debugging = isTRUE(debug);
SEXP keep, loglikelihood, metadata, complete_nodes;

  /* allocate the return value: a vector with length equal to the sample size if
   * we are returning the log-likelihood of each observation, or a vector of
   * length 1 if we are returning the log-likelihood of the whole data set. */
  if (by) {

    PROTECT(loglikelihood = allocVector(REALSXP, ndata));
    loglik = REAL(loglikelihood);
    memset(loglik, '\0', sizeof(double) * ndata);

  }/*THEN*/
  else {

    PROTECT(loglikelihood = ScalarReal(0));
    loglik = Calloc1D(ndata, sizeof(double));

  }/*ELSE*/

  /* find out which nodes to compute the log-likelihood for, zero-indexed. */
  PROTECT(keep = match(keep_nodes, getAttrib(fitted, R_NamesSymbol), 0));

  /* extract the metadata from the data. */
  PROTECT(metadata = getAttrib(data, BN_MetaDataSymbol));
  PROTECT(complete_nodes = getListElement(metadata, "complete.nodes"));

  if (bn.type == DNET || bn.type == ONET || bn.type == DONET) {

    if (debugging)
      Rprintf("> computing the log-likelihood of a discrete network.\n");

    ddata dt = ddata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep);

    if (by) {

      bysample_discrete_loglikelihood(bn, dt, loglik, debugging);

    }/*THEN*/
    else {

      NUM(loglikelihood) =
        data_discrete_loglikelihood(bn, dt, propagate, debugging);

    }/*ELSE*/

    FreeDDT(dt);

  }/*THEN*/
  else if (bn.type == GNET) {

    if (debugging)
      Rprintf("> computing the log-likelihood of a Gaussian network.\n");

    gdata dt = gdata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep);

    if (by) {

      bysample_gaussian_loglikelihood(bn, dt, loglik, debugging);

    }/*THEN*/
    else {

      NUM(loglikelihood) =
        data_gaussian_loglikelihood(bn, dt, loglik, propagate, debugging);

    }/*ELSE*/

    FreeGDT(dt);

  }/*THEN*/
  else if (bn.type == CGNET) {

    if (debugging)
      Rprintf("> computing the log-likelihood of a conditional Gaussian network.\n");

    cgdata dt = cgdata_from_SEXP(data, 0, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete_nodes, keep);

    if (by) {

      bysample_clgaussian_loglikelihood(bn, dt, loglik, debugging);

    }/*THEN*/
    else {

      NUM(loglikelihood) =
        data_clgaussian_loglikelihood(bn, dt, loglik, propagate, debugging);

    }/*ELSE*/

    FreeCGDT(dt);

  }/*THEN*/
  else {

    error("unknown network type, unable to compute the log-likelihood.");

  }/*ELSE*/

  if (!by)
    Free1D(loglik);
  FreeFittedBN(bn);
  UNPROTECT(4);

  return loglikelihood;

}/*DATA_LOGLIKELIHOOD2*/

