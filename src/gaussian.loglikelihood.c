#include "include/rcore.h"
#include "include/blas.h"
#include "include/dataframe.h"
#include "include/globals.h"

double glik(SEXP x, double *nparams) {

int i = 0, num = length(x);
double *xx = REAL(x), xm = 0, res = 0;
long double sd = 0;

  /* compute the mean, */
  for (i = 0; i < num; i++)
    xm += xx[i];
  xm /= num;

  /* compute the standard deviation. */
  for (i = 0; i < num; i++)
    sd += (xx[i] - xm) * (xx[i] - xm);
  sd = sqrt(sd / (num - 1));

  /* compute the log-likelihood (singular models haze zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else
    for (i = 0; i < num; i++)
      res += dnorm(xx[i], xm, (double)sd, TRUE);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = 1;

  return res;

}/*GLIK*/

double cglik(SEXP x, SEXP data, SEXP parents, double *nparams) {

int i = 0, nrow = length(x), ncol = length(parents) + 1;
double *xx = REAL(x), *qr = NULL, res = 0, *fitted = NULL;
long double sd = 0;
SEXP qr_x;

  /* prepare the data for the QR decomposition. */
  PROTECT(qr_x = qr_matrix(data, parents));
  qr = REAL(qr_x);

  /* allocate the fitted values. */
  fitted = Calloc1D(nrow, sizeof(double));

  /* estimate OLS using the QR decomposition. */
  c_qr_ols(qr, xx, nrow, ncol, fitted, &sd);

  /* compute the log-likelihood (singular models haze zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else
    for (i = 0; i < nrow; i++)
      res += dnorm(xx[i], fitted[i], (double)sd, TRUE);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = ncol;

  Free1D(fitted);

  UNPROTECT(1);

  return res;

}/*CGLIK*/

double loglik_gnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel) {

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

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_GNODE*/
