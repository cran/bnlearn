#include "common.h"
#include <Rmath.h>

/* gaussian loss of a node without parents. */
SEXP gloss(SEXP fitted, SEXP data) {

int i = 0, ndata = LENGTH(data);
double *mean = NULL, *sd = NULL, *res = NULL, *x = REAL(data);
SEXP result;

  /* get the coefficient of the linear regression and the standard deviation
   * of the residuals. */
  mean = REAL(getListElement(fitted, "coefficients"));
  sd = REAL(getListElement(fitted, "sd"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* compute the log-likelihood of the data. */
  for (i = 0; i < ndata; i++)
    *res += dnorm(x[i], *mean, *sd, TRUE);

  /* switch to the negentropy. */
  *res /= -ndata;

  UNPROTECT(1);

  return result;

}/*GLOSS*/

/* gaussian loss of a node with one or more parents. */
SEXP cgloss(SEXP fitted, SEXP data)  {

int i = 0, j = 0, ndata = LENGTH(VECTOR_ELT(data, 0)), ncoefs = LENGTH(data);
double mean = 0;
double *res = NULL, *coefs = NULL, *sd = NULL;
double **columns = NULL;
SEXP result;

  /* get the coefficient of the linear regression and the standard deviation
   * of the residuals. */
  coefs = REAL(getListElement(fitted, "coefficients"));
  sd = REAL(getListElement(fitted, "sd"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* dereference the columns of the data frame. */
  columns = (double **) alloc1dpointer(ncoefs);
  for (i = 0; i < ncoefs; i++)
    columns[i] = REAL(VECTOR_ELT(data, i));

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    mean = coefs[0];

    for (j = 1; j < ncoefs; j++)
      mean += columns[j][i] * coefs[j];

    /* compute the log-likelihood of this observation. */
    *res += dnorm(columns[0][i], mean, *sd, TRUE);

  }/*FOR*/

  /* switch to the negentropy. */
  *res /= -ndata;

  UNPROTECT(1);

  return result;

}/*CGLOSS*/

/* multinomial loss of a node without parents. */
SEXP dloss(SEXP fitted, SEXP data, SEXP node) {

int i = 0, ndata = LENGTH(data), dropped = 0;
int *x = INTEGER(data);
double *prob = NULL, *res = NULL;
double logprob = 0;
SEXP result;

  /* get the probabilities of the multinomial distribution. */
  prob = REAL(getListElement(fitted, "prob"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  for (i = 0; i < ndata; i++) {

    logprob = log(prob[x[i] - 1]);

    if (!R_FINITE(logprob) || ISNAN(logprob))
      dropped++;
    else
      *res += logprob;

  }/*FOR*/

  /* switch to the negentropy. */
  *res /= -(ndata - dropped);

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because the corresponding probabilities for node %s were 0 or NaN.", dropped, CHAR(STRING_ELT(node, 0)));

  UNPROTECT(1);

  return result;

}/*DLOSS*/

/* multinomial loss of a node with one or more parents. */
SEXP cdloss(SEXP fitted, SEXP data, SEXP parents, SEXP node) {

int i = 0, ndata = LENGTH(data), nrows = 0, dropped = 0;
int *x = INTEGER(data), *configs = INTEGER(parents);
double *prob = NULL, *res = NULL;
double logprob = 0;
SEXP temp, result;

  /* get the probabilities of the multinomial distribution. */
  temp = getListElement(fitted, "prob");
  nrows = INT(getAttrib(temp, R_DimSymbol));
  prob = REAL(temp);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  for (i = 0; i < ndata; i++) {

    logprob = log(prob[CMC(x[i] - 1, configs[i] - 1, nrows)]);

    if (!R_FINITE(logprob) || ISNAN(logprob))
      dropped++;
    else
      *res += logprob;

  }/*FOR*/

  /* switch to the negentropy. */
  *res /= -(ndata - dropped);

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because the corresponding probabilities for node %s were 0 or NaN.", dropped, CHAR(STRING_ELT(node, 0)));

  UNPROTECT(1);

  return result;

}/*CDLOSS*/

/* classification error of a single node as a loss function. */
SEXP class_err(SEXP reference, SEXP predicted) {

int i = 0, dropped = 0, ndata = LENGTH(reference);
int *r = INTEGER(reference), *p = INTEGER(predicted);
double *res = NULL;
SEXP result;

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* count how many elements differ (this assumes the levels of the two factors
   * are the same and in the same order); NAs are dropped. */
  for (i = 0; i < ndata; i++) {

    if ((r[i] == NA_INTEGER) || (p[i] == NA_INTEGER))
      dropped++;
    else if ((r[i] != p[i]))
      (*res)++;

  }/*FOR*/

  /* rescale into a probability. */
  *res /= (ndata - dropped);

  /* print a warning if data were dropped. */
  if (dropped > 0)
    warning("%d observations were dropped because of missing values.", dropped);

  UNPROTECT(1);

  return result;

}/*CLASS_ERR*/
