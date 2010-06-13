#include "common.h"
#include <R_ext/Applic.h>
#include <Rmath.h>

#define TOTAL_VARIANCE        1
#define GENERALIZED_VARIANCE  2
#define FROBENIUS_NORM        3

static void mvber_var(short int *buffer, double *gen, int k, int *B);
static double _tvar(double *matrix, int *n);
static double _nvar(double *matrix, int *n);

/* test the variance of a multivariate Bernoulli distribution against
 * its worst case value. */
SEXP mvber_variance_test (SEXP var, SEXP replications, SEXP samples,
    SEXP method, SEXP debug) {

int i = 0, r = 0, k = nrows(var);
short int *buffer = NULL;
double *gen = NULL, *res = NULL, tmp = 0;
int *B = INTEGER(samples), *R = INTEGER(replications), *debuglevel = LOGICAL(debug);
SEXP result, genvar;

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, 2));
  res = REAL(result);
  memset(res, '\0', 2 * sizeof(double));

  /* allocate the buffer for the bernoulli random values. */
  buffer = allocstatus(k);

  /* allocate and dereference the generated covariance matrix. */
  PROTECT(genvar = allocMatrix(REALSXP, k, k));
  gen = REAL(genvar);

  /* initialize the random number generator. */
  GetRNGstate();

  switch(INT(method)) {

    case TOTAL_VARIANCE:

      res[0] = _tvar(REAL(var), &k);

      if (*debuglevel > 0)
        Rprintf("* observed value of the test statistic is %lf.\n", res[0]);

      for (r = 0; r < *R; r++) {

        mvber_var(buffer, gen, k, B);

        tmp = _tvar(gen, &k);

        if (*debuglevel > 0)
          Rprintf("  > test statistic is equal to %lf in permutation %d.\n", tmp, r + 1);

        if (res[0] >= tmp)
          res[1] += 1;

      }/*FOR*/

      break;

    case GENERALIZED_VARIANCE:

      /* use the normalized determinant to improve numeric precision;
       * otherwise in many cases the determinant is equal to zero for
       * all the permutations due to the (limited) floating point
       * precision. */
      res[0] = NUM(r_det(var, 4));

      if (*debuglevel > 0)
        Rprintf("* observed value of the test statistic is %lf.\n", res[0]);

      for (r = 0; r < *R; r++) {

        mvber_var(buffer, gen, k, B);

        /* rescale the generated covariance matrix.*/
        for (i = 0; i < k * k; i++)
          gen[i] *= 4;

        tmp = c_det(gen, &k);

        if (*debuglevel > 0)
          Rprintf("  > test statistic is equal to %lf in permutation %d.\n", tmp, r + 1);

        if (res[0] >= tmp)
          res[1] += 1;

      }/*FOR*/

      /* return the real value of the determinant. */
      res[0] = NUM(r_det(var, 1));

      break;

    case FROBENIUS_NORM:

      res[0] = _nvar(REAL(var), &k);

      if (*debuglevel > 0)
        Rprintf("* observed value of the test statistic is %lf.\n", res[0]);

      for (r = 0; r < *R; r++) {

        mvber_var(buffer, gen, k, B);

        tmp = _nvar(gen, &k);

        if (*debuglevel > 0)
          Rprintf("  > test statistic is equal to %lf in permutation %d.\n", tmp, r + 1);

        if (res[0] <= tmp)
          res[1] += 1;

      }/*FOR*/

      break;

  }/*SWITCH*/

  PutRNGstate();

  /* rescale the counter into a significance value. */
  res[1] = res[1] / (double)(*R);

  UNPROTECT(2);

  return result;

}/*MVBER_VARIANCE_TEST*/

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

/* generate a covariance matrix from B multivariate Bernoulli random vectors. */
static void mvber_var(short int *buffer, double *gen, int k, int *B) {

int b = 0, i = 0, j = 0;

  /* clean up the previous covariance matrix. */
  memset(gen, '\0', k * k * sizeof(double));

  for (b = 0; b < *B; b++) {

    /* clean up the buffer. */
    memset(buffer, '\0', k * sizeof(short int));

    /* generate the multivariate Bernoulli random vector. */
    for (i = 0; i < k; i++)
      buffer[i] = unif_rand() > 0.5;

    /* update the co-presence counters in genvar. */
    for (i = 0; i < k; i++)
      for (j = i; j < k; j++)
        if ((buffer[i] == 1) && (buffer[j] == 1)) {

          if (i == j) {

            gen[CMC(i, i, k)]++;

          }/*THEN*/
          else {

            gen[CMC(i, j, k)]++;
            gen[CMC(j, i, k)]++;

          }/*ELSE*/

        }/*THEN*/

  }/*FOR*/

  /* make genvar a real covariance matrix. */
  /* rescale all counts to probabilties. */
  for (i = 0; i < k; i++)
    for (j = i; j < k; j++) {

      gen[CMC(i, j, k)] = gen[CMC(i, j, k)]/(double)(*B);

      if (i != j)
        gen[CMC(j, i, k)] = gen[CMC(i, j, k)];

    }/*FOR*/

  /* compute the covariances first. otherwise the first order probabilities
   * are converted into variances at the same time and mess everything up. */
  for (i = 0; i < k; i++)
    for (j = i+1; j < k; j++)
      gen[CMC(j, i, k)] = gen[CMC(i, j, k)] = gen[CMC(i, j, k)] - gen[CMC(i, i, k)] * gen[CMC(j, j, k)];

  /* compute the variances. */
  for (i = 0; i < k; i++)
    gen[CMC(i, i, k)] = gen[CMC(i, i, k)] - gen[CMC(i, i, k)] * gen[CMC(i, i, k)];

}/*MVBER_VAR*/

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

      res[0] = _nvar(REAL(var), &k);
      res[1] = NA_REAL;
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

static double _nvar(double *matrix, int *n) {

int i = 0, j = 0;
double res = 0;

  for (i = 0; i < *n; i++)
    for (j = i+1; j < *n; j++)
      res += 2 * matrix[CMC(i, j, *n)] * matrix[CMC(i, j, *n)];

    for (i = 0; i < *n; i++)
      res += (matrix[CMC(i, i, *n)] - 0.25) * (matrix[CMC(i, i, *n)] - 0.25);

  return res;

}/*_NVAR*/

