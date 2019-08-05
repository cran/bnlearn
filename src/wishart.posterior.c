#include "include/rcore.h"
#include "include/blas.h"
#include "include/covariance.h"
#include "include/data.frame.h"
#include "include/data.table.h"
#include "include/matrix.h"
#include "include/scores.h"

double wpost(double *xx, int data_ncols, int nobs, double alpha_mu, double *nu,
    double alpha_w) {

double r = 0, t = 0, xm = 0, xsse = 0;
long double logprob = 0;

  /* first term. */
  logprob = 0.5 * (log(alpha_mu) - log(nobs + alpha_mu));

  /* Gamma_l ratio in the second term. */
  logprob += lgammafn(0.5 * (nobs + alpha_w - data_ncols + 1)) -
             lgammafn(0.5 * (alpha_w - data_ncols+ 1));

  /* leftover from the second term. */
  logprob -= 0.5 * nobs * log(M_PI);

  /* third term, numerator. */
  t = alpha_mu * (alpha_w - data_ncols - 1) / (alpha_mu + 1);
  logprob += 0.5 * (alpha_w - data_ncols + 1) * log(t);

  /* third term, denominator. */
  xm = c_mean(xx, nobs);
  xsse = c_sse(xx, xm, nobs);

  r = t + xsse + (nobs * alpha_mu) / (nobs + alpha_mu) * (xm - *nu) * (xm - *nu);

  logprob -= 0.5 * (nobs + alpha_w - data_ncols + 1) * log(r);

  return (double)logprob;

}/*WPOST*/

double cwpost(double *xx, SEXP z, int data_ncols, int nobs, double alpha_mu,
    double *nu, double alpha_w) {

int i = 0, j = 0, p = length(z);
double t = 0;
long double logprob = 0;

  /* first term. */
  logprob = 0.5 * (log(alpha_mu) - log(nobs + alpha_mu));

  /* Gamma_l ratio in the second term. */
  logprob += lgammafn(0.5 * (nobs + alpha_w - data_ncols + p + 1)) -
             lgammafn(0.5 * (alpha_w - data_ncols + p + 1));

  /* leftover from the second term. */
  logprob -= 0.5 * nobs * log(M_PI);

  /* third term, ratio of the determinants of the prior T matrices. */
  t = alpha_mu * (alpha_w - data_ncols - 1) / (alpha_mu + 1);

  logprob += 0.5 * (alpha_w - data_ncols + p + 1) * (p + 1) * log(t) -
             0.5 * (alpha_w - data_ncols + p) * (p) * log(t);

  /* third term, ratio of the determinants of the posterior R matrices. */
  gdata dtx = gdata_from_SEXP(z, 1);
  dtx.col[0] = xx;
  gdata_cache_means(&dtx, 0);
  covariance R = new_covariance(dtx.m.ncols, FALSE);
  covariance Rtilde = new_covariance(dtx.m.ncols - 1, FALSE);

  /* compute the rescaled covariance matrix. */
  c_covmat(dtx.col, dtx.mean, dtx.m.nobs, dtx.m.ncols, R, 0);
  for (i = 0; i < R.dim * R.dim; i++)
    R.mat[i] *= nobs - 1;

  /* add the prior matrix, diagonal with elements t. */
  for (i = 0; i < R.dim; i++)
    R.mat[CMC(i, i, R.dim)] += t;

  /* add the outer product of the difference in observed and prior means. */
  for (i = 0; i < R.dim; i++)
    for (j = 0; j < R.dim; j++) {

      R.mat[CMC(i, j, R.dim)] += (nobs * alpha_mu) / (nobs + alpha_mu) *
                                  (dtx.mean[i] - nu[i]) * (dtx.mean[j] - nu[j]);

  }/*FOR*/

  /* subset the rescaled covariance matrix to include only the parents. */
  covariance_drop_variable(&R, &Rtilde, 0);

  logprob += 0.5 * (nobs + alpha_w - data_ncols + p) *
               log(c_det(Rtilde.mat, Rtilde.dim));
  logprob -= 0.5 * (nobs + alpha_w - data_ncols + p + 1) *
               log(c_det(R.mat, R.dim));

  FreeGDT(dtx, FALSE);
  FreeCOV(R);
  FreeCOV(Rtilde);

  return (double)logprob;

}/*CWPOST*/

double wishart_node(SEXP target, SEXP x, SEXP data, SEXP iss_mu, SEXP nu,
    SEXP iss_w, SEXP prior, SEXP beta, bool debugging) {

int nobs = length(VECTOR_ELT(data, 0)), ncols = length(data);
char *t = (char *)CHAR(STRING_ELT(target, 0));
double prob = 0, prior_prob = 0;
SEXP nodes, node_t, parents, parent_vars, data_t, nu_subset;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
  /* compute the prior probability component for the node. */
  prior_prob = graph_prior_prob(prior, target, beta, nodes, debugging);
  /* extract the elements of the prior mean vector. */
  PROTECT(nu_subset = subset_by_name(nu, 2, target, parents));

  if (length(parents) == 0) {

    prob = wpost(REAL(data_t), ncols, nobs, NUM(iss_mu), REAL(nu_subset),
             NUM(iss_w));

  }/*THEN*/
  else {

    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    /* compute the marginal likelihood. */
    prob = cwpost(REAL(data_t), parent_vars, ncols, nobs, NUM(iss_mu),
             REAL(nu_subset), NUM(iss_w));

    UNPROTECT(1);

  }/*ELSE*/

  if (debugging) {

    Rprintf("  > (log)prior probability is %lf.\n", prior_prob);
    Rprintf("  > (log)posterior density is %lf.\n", prob);

  }/*THEN*/

  /* add the (log)prior to the marginal (log)likelihood to get the (log)posterior. */
  prob += prior_prob;

  UNPROTECT(2);

  return prob;

}/*WISHART_NODE*/
