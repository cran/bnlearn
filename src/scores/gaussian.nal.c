#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/data.frame.h"
#include "../minimal/common.h"
#include "../include/globals.h"
#include "../core/moments.h"
#include "../math/linear.algebra.h"

double glik_incomplete(SEXP x, double k) {

int i = 0, num = length(x), nc = 0;
double *xx = REAL(x), res  = 0, xm = 0, sd = 0;

  c_ols(NULL, xx, num, 0, NULL, NULL, &xm, &sd, &nc, TRUE);

  /* compute the log-likelihood. */
  if ((sd < MACHINE_TOL) || (nc == 0)) {

    /* singular distributions and those with no locally-complete observations
     * are scored to prevent their selection in structure learning. */
    res = R_NegInf;

  }/*THEN*/
  else {

    for (i = 0; i < num; i++)
      if (!ISNAN(xx[i]))
        res += dnorm(xx[i], xm, (double)sd, TRUE);

    /* average over the locally-complete data... */
    res /= nc;
    /* and add the penalty term, scaled by the original sample size. */
    res -= k / num * 2;

  }/*ELSE*/

  return res;

}/*GLIK_INCOMPLETE*/

double cglik_incomplete(SEXP x, SEXP data, SEXP parents, double k) {

int i = 0, nrow = length(x), ncol = length(parents), nc = 0;
double *xx = REAL(x), **dd = NULL, res = 0, *fitted = NULL, sd = 0;
SEXP data_x;

  /* dereference the data. */
  PROTECT(data_x = c_dataframe_column(data, parents, FALSE, FALSE));
  dd = Calloc1D(ncol, sizeof(double *));
  for (i = 0; i < ncol; i++)
    dd[i] = REAL(VECTOR_ELT(data_x, i));
  /* allocate the fitted values. */
  fitted = Calloc1D(nrow, sizeof(double));

  c_ols(dd, xx, nrow, ncol, fitted, NULL, NULL, &sd, &nc, TRUE);

  /* compute the log-likelihood. */
  if ((sd < MACHINE_TOL) || (nc == 0)) {

    /* singular distributions and those with no locally-complete observations
     * are scored to prevent their selection in structure learning. */
    res = R_NegInf;

  }/*THEN*/
  else {

    for (i = 0; i < nrow; i++)
      if (!ISNAN(fitted[i]) && !ISNAN(xx[i]))
        res += dnorm(xx[i], fitted[i], (double)sd, TRUE);

    /* average over the locally-complete data... */
    res /= nc;
    /* and add the penalty term, scaled by the original sample size. */
    res -= k / nrow * (ncol + 2);

  }/*ELSE*/

  Free1D(fitted);
  Free1D(dd);

  UNPROTECT(1);

  return res;

}/*CGLIK_INCOMPLETE*/

double nal_gnode(SEXP target, SEXP x, SEXP data, double k, bool debugging) {

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
    loglik = glik_incomplete(data_t, k);
  else
    loglik = cglik_incomplete(data_t, data, parents, k);

  if (debugging)
    Rprintf("  > log-likelihood is %lf.\n", loglik);

  return loglik;

}/*NAL_GNODE*/
