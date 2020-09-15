#include "include/rcore.h"
#include "include/blas.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"

double pgnode(SEXP x, SEXP new_x, double *nparams) {

int i = 0, num = length(x), num2 = length(new_x);
double *xx = REAL(x), *xx2 = REAL(new_x), xm = 0, res = 0;
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

  /* compute the log-likelihood (singular models have zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else
    for (i = 0; i < num2; i++)
      res += dnorm(xx2[i], xm, (double)sd, TRUE);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = 1;

  return res;

}/*PGNODE*/

double cpgnode(SEXP x, SEXP x2, SEXP data, SEXP newdata, SEXP parents,
    double *nparams) {

int i = 0, j = 0, nrow = length(x), nrow2 = length(x2), ncol = length(parents);
double *xx = REAL(x), *xx2 = REAL(x2), **dd = NULL, **dd2 = NULL;
double res = 0, sd = 0, *beta = NULL, fitted = 0;
SEXP data_x, new_x;

  /* dereference the data. */
  PROTECT(data_x = c_dataframe_column(data, parents, FALSE, FALSE));
  dd = Calloc1D(ncol, sizeof(double *));
  for (i = 0; i < ncol; i++)
    dd[i] = REAL(VECTOR_ELT(data_x, i));

  PROTECT(new_x = c_dataframe_column(newdata, parents, FALSE, FALSE));
  dd2 = Calloc1D(ncol, sizeof(double *));
  for (i = 0; i < ncol; i++)
    dd2[i] = REAL(VECTOR_ELT(new_x, i));

  beta = Calloc1D(ncol + 1, sizeof(double));

  c_ols(dd, xx, nrow, ncol, NULL, NULL, beta, &sd, FALSE);

  /* this estimate is not unbiased; the denominator is set to match the
   * corresponding denominator in the root nodes to ensure score equivalence. */
  sd = sd * sqrt((double) (nrow - ncol - 1) / (double) (nrow - 1));

  /* compute the log-likelihood (singular models have zero density). */
  if (sd < MACHINE_TOL)
    res = R_NegInf;
  else {

    for (i = 0; i < nrow2; i++) {

      for (j = 1, fitted = beta[0]; j < ncol + 1; j++)
        fitted += beta[j] * dd2[j - 1][i];

      res += dnorm(xx2[i], fitted, sd, TRUE);

    }/*FOR*/

  }/*ELSE*/

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = ncol + 1;

  Free1D(beta);
  Free1D(dd);
  Free1D(dd2);

  UNPROTECT(2);

  return res;

}/*CPGNODE*/

double predictive_loglik_gnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, int debugging) {

double loglik = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, new_t;

  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
  PROTECT(new_t = c_dataframe_column(newdata, target, TRUE, FALSE));

  /* compute the log-likelihood. */
  if (length(parents) == 0)
    loglik = pgnode(data_t, new_t, nparams);
  else
    loglik = cpgnode(data_t, new_t, data, newdata, parents, nparams);

  if (debugging > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  UNPROTECT(2);

  return loglik;

}/*PREDICTIVE_LOGLIK_GNODE*/
