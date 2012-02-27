#include "common.h"
#include <R_ext/Applic.h>
#include <Rmath.h>

#define TOTAL_VARIANCE        1
#define GENERALIZED_VARIANCE  2
#define FROBENIUS_NORM        3
#define FROBENIUS_NORM_K      4

static double _tvar(double *matrix, int *n);
static double _nvar(double *matrix, int *n, double lambda);

/* compute the joint probabilities of each pair of arcs. */
SEXP mvber_joint_counters(SEXP arcs, SEXP nodes, SEXP joint, SEXP debug) {

int i = 0, j = 0, k = 0, *coords = NULL;
int nnodes = LENGTH(nodes), narcs = LENGTH(arcs)/2, jsize = nrows(joint);
int bufsize = UPTRI3_MATRIX(nnodes);
short int *buffer = NULL;
double *p = REAL(joint);
SEXP try;

  /* allocate the buffer for the arc identifiers. */
  buffer = allocstatus(bufsize);

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  /* use the UPTRI3() macro as a unique identifier of each arc (modulo its
   * direction) and store the resulting hashes in the buffer. */
  for (k = 0; k < narcs; k++)
    buffer[UPTRI3(coords[k], coords[k + narcs], nnodes)] = 1;

  /* print arcs' mappings, which are essential to get the moments right. */
  if (isTRUE(debug)) {

    for (k = 0; k < narcs; k++) {

      Rprintf("  > learned arc (%s, %s) -> mapped to (%d, %d) = %d\n",
        CHAR(STRING_ELT(arcs, k)),
        CHAR(STRING_ELT(arcs, k + narcs)),
        coords[k], coords[k + narcs],
        UPTRI3(coords[k], coords[k + narcs], nnodes));

    }/*FOR*/

  }/*THEN*/

  /* unprotect coords, it's no longer needed. */
  UNPROTECT(1);

  /* increment the joint counters in the joint matrix. */
  for (i = 0; i < bufsize; i++) {

    for (j = i; j < bufsize; j++)  {

      if ((buffer[i] == 1) && (buffer[j] == 1)) {

        if (i == j) {

          p[CMC(i, i, jsize)]++;

        }/*THEN*/
        else {

          p[CMC(i, j, jsize)]++;
          p[CMC(j, i, jsize)]++;

        }/*ELSE*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return joint;

}/*MVBER_JOINT_COUNTERS*/

/* convert the joint probabilities to the first and second moments of
 * the multivariate Bernoulli distribution. */
SEXP mvber_moments(SEXP joint, SEXP R, SEXP labels) {

int i = 0, j = 0, dim = nrows(joint);
double *p = NULL, *n = REAL(joint), *r = REAL(R);
SEXP res, names, first, dimnames;

  /* allocate the return value. */
  PROTECT(res = allocVector(VECSXP, 2));

  /* allocate and set the elements' names. */
  PROTECT(names = allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("expected"));
  SET_STRING_ELT(names, 1, mkChar("covariance"));
  setAttrib(res, R_NamesSymbol, names);

  /* rescale all counts to probabilties. */
  for (i = 0; i < dim; i++)
    for (j = i; j < dim; j++) {

      n[CMC(i, j, dim)] = n[CMC(i, j, dim)] / (*r);

      if (i != j)
        n[CMC(j, i, dim)] = n[CMC(i, j, dim)];

    }/*FOR*/

  /* allocate and dereference the expected value. */
  PROTECT(first = allocVector(REALSXP, dim));
  p = REAL(first);

  /* initialize the expected value. */
  for (i = 0; i < dim; i++)
    p[i] = n[CMC(i, i, dim)];

  /* set the labels on the probability vector. */
  setAttrib(first, R_NamesSymbol, labels);

  /* make joint a true covariance matrix. */

  /* compute the covariances first. otherwise the first order probabilities
   * are converted into variances at the same time and mess everything up. */
  for (i = 0; i < dim; i++)
    for (j = i+1; j < dim; j++)
      n[CMC(j, i, dim)] = n[CMC(i, j, dim)] = n[CMC(i, j, dim)] - n[CMC(i, i, dim)] * n[CMC(j, j, dim)];

  /* compute the variances. */
  for (i = 0; i < dim; i++)
    n[CMC(i, i, dim)] = n[CMC(i, i, dim)] - n[CMC(i, i, dim)] * n[CMC(i, i, dim)];

  /* set the labels on the covariance matrix. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, labels);
  SET_VECTOR_ELT(dimnames, 1, labels);
  setAttrib(joint, R_DimNamesSymbol, dimnames);

  /* store the expected value and the covariance matrix in the return value. */
  SET_VECTOR_ELT(res, 0, first);
  SET_VECTOR_ELT(res, 1, joint);

  UNPROTECT(4);

  return res;

}/*MVBER_MOMENTS*/

/* descriptive statistics of network's variability. */
SEXP mvber_variance(SEXP var, SEXP method) {

int k = nrows(var);
double *res = NULL;
SEXP result;


  /* allocate the return value. */
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);

  switch(INT(method)) {

    case TOTAL_VARIANCE:

      res[0] = _tvar(REAL(var), &k);
      res[1] = 4 * res[0] / k;
      break;

    case GENERALIZED_VARIANCE:

      res[0] = NUM(r_det(var, 1));
      res[1] = NUM(r_det(var, 4));
      break;

    case FROBENIUS_NORM:

      res[0] = _nvar(REAL(var), &k, 0.25);
      res[1] = NA_REAL;
      break;

    case FROBENIUS_NORM_K:

      res[0] = _nvar(REAL(var), &k, 0.25 * k);
      res[1] = (res[0] - 0.0625 * k * (k - 1) * (k - 1)) / (0.0625 * k * k * k - 0.0625 * k * (k - 1) * (k - 1));
      break;

    default:

      error("unknown decriptive statistic of the network's structure variability");
      break;

  }/*SWITCH*/

  UNPROTECT(1);

  return result;

}/*MVBER_VARIANCE*/

static double _tvar(double *matrix, int *n) {

int i = 0;
double res = 0;

  for (i = 0; i < *n; i++)
    res += matrix[CMC(i, i, *n)];

  return res;

}/*_TVAR*/

static double _nvar(double *matrix, int *n, double lambda) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < *n; i++)
    for (j = i+1; j < *n; j++)
      res += 2 * matrix[CMC(i, j, *n)] * matrix[CMC(i, j, *n)];

    for (i = 0; i < *n; i++)
      res += (matrix[CMC(i, i, *n)] - lambda) * (matrix[CMC(i, i, *n)] - lambda);

  return res;

}/*_NVAR*/

