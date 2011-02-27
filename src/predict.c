#include "common.h"

/* predict the value of a gaussian node without parents. */
SEXP gpred(SEXP fitted, SEXP data, SEXP debug) {

int i = 0, *ndata = INTEGER(data), *debuglevel = LOGICAL(debug);
double *mean = NULL, *res = NULL;
SEXP result;

  /* get the (only) coefficient of the linear regression. */
  mean = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, *ndata));
  res = REAL(result);

  /* copy the mean in the return value. */
  for (i = 0; i < *ndata; i++)
    res[i] = *mean;

  if (*debuglevel > 0) {

    Rprintf("  > prediction for all observations is %lf.\n", *mean);

  }/*THEN*/

  UNPROTECT(1);

  return result;

}/*GPRED*/

/* predict the value of a gaussian node with one or more parents. */
SEXP cgpred(SEXP fitted, SEXP data, SEXP debug)  {

int i = 0, j = 0, ndata = LENGTH(VECTOR_ELT(data, 0)), ncols = LENGTH(data);
int *debuglevel = LOGICAL(debug);
double *res = NULL, *coefs = NULL;
double **columns = NULL;
SEXP result;

  /* get the coefficient of the linear regression. */
  coefs = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, ndata));
  res = REAL(result);

  /* dereference the columns of the data frame. */
  columns = Calloc(ncols, double *);
  for (i = 0; i < ncols; i++)
    columns[i] = REAL(VECTOR_ELT(data, i));

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    res[i] = coefs[0];

    for (j = 0; j < ncols; j++)
      res[i] += columns[j][i] * coefs[j + 1];

    if (*debuglevel > 0) {

      Rprintf("  > prediction for observation %d is %lf with predictor:\n",
        i + 1, res[i]);

      Rprintf("    (%lf) + (%lf) * (%lf)", coefs[0], columns[0][i], coefs[1]);
      for (j = 1; j < ncols; j++) 
        Rprintf(" + (%lf) * (%lf)", columns[j][i], coefs[j + 1]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  Free(columns);

  UNPROTECT(1);

  return result;

}/*CGPRED*/

/* predict the value of a discrete node without parents. */
SEXP dpred(SEXP fitted, SEXP data, SEXP debug) {

int i = 0, imax = 0, ndata = LENGTH(data);
int *res = NULL, *debuglevel = LOGICAL(debug);
double *prob = NULL; 
SEXP ptab, result, tr_levels = getAttrib(data, R_LevelsSymbol);

  /* get the probabilities of the multinomial distribution. */
  ptab = getListElement(fitted, "prob");
  prob = REAL(ptab);

  /* find out the mode. */
  imax = which_max(prob, LENGTH(ptab));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  /* copy the index of the mode in the return value. */
  for (i = 0; i < ndata; i++)
    res[i] = imax;

  if (*debuglevel > 0) {

    if (res[0] == NA_INTEGER)
      Rprintf("  > prediction for all observations is NA with probabilities:\n");
    else
      Rprintf("  > prediction for all observations is %s with probabilities:\n",
      CHAR(STRING_ELT(tr_levels, res[0] - 1)));

    Rprintf("  ");
    for (i = 0; i < LENGTH(ptab); i++)
      Rprintf("  %lf", prob[i]);
    Rprintf("\n");

  }/*THEN*/

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*DPRED*/

/* predict the value of a discrete node with one or more parents. */
SEXP cdpred(SEXP fitted, SEXP data, SEXP parents, SEXP debug) {

int i = 0, ndata = LENGTH(data), nrows = 0, ncols = 0;
int *configs = INTEGER(parents), *mode = NULL, *res = NULL;
int *debuglevel = LOGICAL(debug);
double *prob = NULL;
SEXP temp, result, tr_levels = getAttrib(data, R_LevelsSymbol);

  /* get the probabilities of the multinomial distribution. */
  temp = getListElement(fitted, "prob");
  nrows = INT(getAttrib(temp, R_DimSymbol));
  ncols = LENGTH(temp) / nrows;
  prob = REAL(temp);

  /* get the mode for each configuration. */
  mode = alloc1dcont(ncols);
  for (i = 0; i < ncols; i++)
    mode[i] = which_max(prob + nrows * i, nrows);

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(INTSXP, ndata));
  res = INTEGER(result);

  /* copy the index of the mode in the return value. */
  for (i = 0; i < ndata; i++) {

    res[i] = mode[configs[i] - 1];

    if (*debuglevel > 0) {

      if (res[i] == NA_INTEGER)
        Rprintf("  > prediction for observation %d is NA with probabilities:\n");
      else
        Rprintf("  > prediction for observation %d is %s with probabilities:\n",
          i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

      Rprintf("  ");
      for (int k = 0; k < nrows; k++) 
        Rprintf("  %lf", (prob + nrows * (configs[i] - 1))[k]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  /* copy the labels and the class from the input data. */
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, getAttrib(data, R_ClassSymbol));

  UNPROTECT(1);

  return result;

}/*CDPRED*/

/* predict the value of the training variable in a naive Bayes classifier. */
SEXP naivepred(SEXP fitted, SEXP data, SEXP training, SEXP prior, SEXP debug) {

int i = 0, j = 0, k = 0, n = 0, nvars = LENGTH(fitted), tr_nlevels = 0;
int *res = NULL, *cpt_nrows = NULL, **ex = NULL, *tr_id = INTEGER(training);
int *debuglevel = LOGICAL(debug);
double **cpt = NULL, *pr = NULL, *scratch = NULL;
SEXP class, temp, tr, tr_levels, result;

  /* cache the pointers to all the variables. */
  ex = (int **) alloc1dpointer(nvars);

  for (i = 0; i < nvars; i++) {

    ex[i] = INTEGER(VECTOR_ELT(data, i));

  }/*FOR*/

  /* get the training variable and its levels. */
  n = LENGTH(VECTOR_ELT(data, 0));
  tr = getListElement(VECTOR_ELT(fitted, *tr_id - 1), "prob");
  tr_levels = VECTOR_ELT(getAttrib(tr, R_DimNamesSymbol), 0);
  tr_nlevels = LENGTH(tr_levels);
  /* get the prior distribution. */
  pr = REAL(prior);
  /* allocate the scratch space used to compute posterior probabilities. */
  scratch = alloc1dreal(tr_nlevels);

  /* cache the pointers to the conditional probability tables. */
  cpt = (double **) alloc1dpointer(nvars);
  cpt_nrows = alloc1dcont(nvars);

  for (i = 0; i < nvars; i++) {

    temp = VECTOR_ELT(fitted, i);
    temp = getListElement(temp, "prob");

    cpt[i] = REAL(temp);
    cpt_nrows[i] = LENGTH(temp) / tr_nlevels;

  }/*FOR*/

  /* allocate the return value. */
  PROTECT(result = allocVector(INTSXP, n));
  res = INTEGER(result);

  /* for each observation... */
  for (i = 0; i < n; i++) {

    /* ... reset the scratch space... */
    for (k = 0; k < tr_nlevels; k++)
      scratch[k] = log(pr[k]);

    /* ... and for each conditional probability table... */
    for (j = 0; j < nvars; j++) {

      /* ... skip the training variable.... */
      if (*tr_id == j + 1)
        continue;

      /* ... and for each row of the conditional probability table... */
      for (k = 0; k < tr_nlevels; k++) {

        /* ... update the posterior probability. */
        scratch[k] += log(cpt[j][CMC(ex[j][i] - 1, k, cpt_nrows[j])]);

      }/*FOR*/

    }/*FOR*/

    res[i] = which_max(scratch, tr_nlevels);

    if (*debuglevel > 0) {

      Rprintf("  > prediction for observation %d is %s with (log-)posterior:\n",
        i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

      Rprintf("  ");
      for (k = 0; k < tr_nlevels; k++) 
        Rprintf("  %lf", scratch[k]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  /* add back the attributes and the class to the return value. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, class);

  UNPROTECT(2);

  return result;

}/*NAIVEPRED*/
