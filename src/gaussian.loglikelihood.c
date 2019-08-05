#include "include/rcore.h"
#include "include/blas.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"

double glik(SEXP x, double *nparams) {

int i = 0, num = length(x);
double *xx = REAL(x), xm = 0, res = 0;
long double sd = 0;

  /* compute the mean, */
  for (i = 0; i < num; i++)
    xm += xx[i];
  xm /= num;

  /* compute the standard deviation. */
  if (num == 0)
    sd = R_NaN;
  else if (num == 1)
    sd = 0;
  else {

    for (i = 0; i < num; i++)
      sd += (xx[i] - xm) * (xx[i] - xm);
    sd = sqrt(sd / (num - 1));

  }/*ELSE*/

  /* compute the log-likelihood (singular models haze zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else
    for (i = 0; i < num; i++)
      res += dnorm(xx[i], xm, (double)sd, TRUE);

  /* we may want to store the number of parameters (one mean, one standard
   * deviation). */
  if (nparams)
    *nparams = 2;

  return res;

}/*GLIK*/

double cglik(SEXP x, SEXP data, SEXP parents, double *nparams) {

int i = 0, nrow = length(x), ncol = length(parents);
double *xx = REAL(x), **dd = NULL, res = 0, *fitted = NULL, sd = 0;
SEXP data_x;

  /* dereference the data. */
  PROTECT(data_x = c_dataframe_column(data, parents, FALSE, FALSE));
  dd = Calloc1D(ncol, sizeof(double *));
  for (i = 0; i < ncol; i++)
    dd[i] = REAL(VECTOR_ELT(data_x, i));
  /* allocate the fitted values. */
  fitted = Calloc1D(nrow, sizeof(double));

  c_ols(dd, xx, nrow, ncol, fitted, NULL, NULL, &sd, FALSE);

  /* compute the log-likelihood (singular models have zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else
    for (i = 0; i < nrow; i++)
      res += dnorm(xx[i], fitted[i], (double)sd, TRUE);

  /* we may want to store the number of parameters (one intercept, one
   * regression coefficient per parent, one residuals standard deviation). */
  if (nparams)
    *nparams = ncol + 2;

  Free1D(fitted);
  Free1D(dd);

  UNPROTECT(1);

  return res;

}/*CGLIK*/

double loglik_gnode(SEXP target, SEXP x, SEXP data, double *nparams,
    bool debugging) {

double loglik = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  data_t = c_dataframe_column(data, target, TRUE, FALSE);

  /* compute the log-likelihood. */
  if (length(parents) == 0)
    loglik = glik(data_t, nparams);
  else
    loglik = cglik(data_t, data, parents, nparams);

  if (debugging)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_GNODE*/
