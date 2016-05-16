#include "include/rcore.h"
#include "include/globals.h"
#include "include/dataframe.h"
#include "include/sampling.h"
#include "include/matrix.h"

/* predict the value of a gaussian node without parents. */
SEXP gpred(SEXP fitted, SEXP ndata, SEXP debug) {

int i = 0, *n = INTEGER(ndata), debuglevel = isTRUE(debug);
double *mean = NULL, *res = NULL;
SEXP result;

  /* get the (only) coefficient of the linear regression. */
  mean = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, *n));
  res = REAL(result);

  /* copy the mean in the return value. */
  for (i = 0; i < *n; i++)
    res[i] = *mean;

  if (debuglevel > 0) {

    Rprintf("  > prediction for all observations is %lf.\n", *mean);

  }/*THEN*/

  UNPROTECT(1);

  return result;

}/*GPRED*/

/* predict the value of a gaussian node with one or more parents. */
SEXP cgpred(SEXP fitted, SEXP parents, SEXP debug)  {

int i = 0, j = 0, ndata = length(VECTOR_ELT(parents, 0)), ncols = length(parents);
int debuglevel = isTRUE(debug);
double *res = NULL, *coefs = NULL;
double **columns = NULL;
SEXP result;

  /* get the coefficient of the linear regression. */
  coefs = REAL(getListElement(fitted, "coefficients"));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, ndata));
  res = REAL(result);

  /* dereference the columns of the data frame. */
  columns = (double **) Calloc1D(ncols, sizeof(double));
  for (i = 0; i < ncols; i++)
    columns[i] = REAL(VECTOR_ELT(parents, i));

  for (i = 0; i < ndata; i++) {

    /* compute the mean value for this observation. */
    res[i] = coefs[0];

    for (j = 0; j < ncols; j++)
      res[i] += columns[j][i] * coefs[j + 1];

    if (debuglevel > 0) {

      Rprintf("  > prediction for observation %d is %lf with predictor:\n",
        i + 1, res[i]);

      Rprintf("    (%lf) + (%lf) * (%lf)", coefs[0], columns[0][i], coefs[1]);
      for (j = 1; j < ncols; j++)
        Rprintf(" + (%lf) * (%lf)", columns[j][i], coefs[j + 1]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  Free1D(columns);

  return result;

}/*CGPRED*/

/* predict the value of a discrete node without parents. */
SEXP dpred(SEXP fitted, SEXP ndata, SEXP prob, SEXP debug) {

int i = 0, nmax = 0, n = INT(ndata), length = 0;
int *res = NULL, debuglevel = isTRUE(debug), include_prob = isTRUE(prob);
int *iscratch = NULL, *maxima = NULL, tr_nlevels = 0;
double *cpt = NULL, *dscratch = NULL, *buf = NULL, *pt = NULL;
SEXP ptab, result, tr_levels, probtab = R_NilValue;

  /* get the probabilities of the multinomial distribution. */
  ptab = getListElement(fitted, "prob");
  length = length(ptab);
  cpt = REAL(ptab);

  /* create the vector of indexes. */
  iscratch = Calloc1D(length, sizeof(int));
  for (i = 0; i < length; i++)
    iscratch[i] = i + 1;

  /* create a scratch copy of the array. */
  buf = Calloc1D(length, sizeof(double));
  dscratch = Calloc1D(length, sizeof(double));
  memcpy(dscratch, cpt, length * sizeof(double));

  /* allocate the array for the indexes of the maxima. */
  maxima = Calloc1D(length, sizeof(int));

  /* find out the mode(s). */
  nmax = all_max(dscratch, length, maxima, iscratch, buf);

  /* allocate and initialize the return value. */
  PROTECT(result = node2df(fitted, n));
  res = INTEGER(result);
  /* copy the levels for use in debuggging. */
  tr_levels = getAttrib(result, R_LevelsSymbol);
  tr_nlevels = length(tr_levels);

  /* allocate and initialize the table of the prediction probabilities. */
  if (include_prob > 0) {

    PROTECT(probtab = allocMatrix(REALSXP, tr_nlevels, n));
    pt = REAL(probtab);
    for (i = 0; i < n; i++)
      memcpy(pt + i * tr_nlevels, cpt, tr_nlevels * sizeof(double));

  }/*THEN*/

  if (nmax == 1) {

    /* copy the index of the mode in the return value. */
    for (i = 0; i < n; i++)
      res[i] = maxima[0];

    if (debuglevel > 0) {

      if (res[0] == NA_INTEGER)
        Rprintf("  > prediction for all observations is NA with probabilities:\n");
      else
        Rprintf("  > prediction for all observations is '%s' with probabilities:\n",
          CHAR(STRING_ELT(tr_levels, res[0] - 1)));

      Rprintf("  ");
      for (i = 0; i < length(ptab); i++)
        Rprintf("  %lf", cpt[i]);
      Rprintf("\n");

    }/*THEN*/

  }/*THEN*/
  else {

    /* break ties: sample with replacement from all the maxima. */
    GetRNGstate();
    SampleReplace(n, nmax, res, maxima);
    PutRNGstate();

    if (debuglevel > 0) {

      Rprintf("  > there are %d levels tied for prediction, applying tie breaking.\n", nmax);
      Rprintf("  > tied levels are:");
      for (i = 0; i < nmax; i++)
        Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[i] - 1)));
      Rprintf(".\n");

    }/*THEN*/

  }/*ELSE*/

  if (include_prob > 0) {

    /* set the levels of the target variable as rownames. */
    setDimNames(probtab, tr_levels, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(2);

  }/*THEN*/
  else {

    UNPROTECT(1);

  }/*ELSE*/

  Free1D(iscratch);
  Free1D(buf);
  Free1D(dscratch);
  Free1D(maxima);

  return result;

}/*DPRED*/

/* predict the value of a discrete node with one or more parents. */
SEXP cdpred(SEXP fitted, SEXP parents, SEXP prob, SEXP debug) {

int i = 0, k = 0, n = length(parents), nrows = 0, ncols = 0;
int *configs = INTEGER(parents), debuglevel = isTRUE(debug);
int tr_nlevels = 0, include_prob = isTRUE(prob);
int *iscratch = NULL, *maxima = NULL, *nmax = NULL, *res = NULL;
double *cpt = NULL, *dscratch = NULL, *buf = NULL, *pt = NULL;
SEXP temp, result, tr_levels, probtab = R_NilValue;

  /* get the probabilities of the multinomial distribution. */
  temp = getListElement(fitted, "prob");
  nrows = INT(getAttrib(temp, R_DimSymbol));
  ncols = length(temp) / nrows;
  cpt = REAL(temp);

  /* create the vector of indexes. */
  iscratch = Calloc1D(nrows, sizeof(int));

  /* create a scratch copy of the array. */
  buf = Calloc1D(nrows, sizeof(double));
  dscratch = Calloc1D(nrows * ncols, sizeof(double));
  memcpy(dscratch, cpt, nrows * ncols * sizeof(double));

  /* allocate the array for the indexes of the maxima. */
  maxima = Calloc1D(nrows * ncols, sizeof(int));

  /* allocate the maxima counters. */
  nmax = Calloc1D(ncols, sizeof(int));

  /* get the mode for each configuration. */
  for (i = 0; i < ncols; i++) {

    /* initialize the vector of indexes. */
    for (k = 0; k < nrows; k++)
      iscratch[k] = k + 1;

    /* find out the mode(s). */
    nmax[i] = all_max(dscratch + i * nrows, nrows, maxima + i * nrows,
                iscratch, buf);

  }/*FOR*/

  /* allocate and initialize the return value. */
  PROTECT(result = node2df(fitted, n));
  res = INTEGER(result);
  /* copy the levels for use in debuggging. */
  tr_levels = getAttrib(result, R_LevelsSymbol);
  tr_nlevels = length(tr_levels);

  /* allocate and initialize the table of the prediction probabilities. */
  if (include_prob > 0) {

    PROTECT(probtab = allocMatrix(REALSXP, tr_nlevels, n));
    pt = REAL(probtab);

  }/*THEN*/

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  /* copy the index of the mode in the return value. */
  for (i = 0; i < n; i++) {

    if (nmax[configs[i] - 1] == 1) {

      res[i] = maxima[CMC(0, configs[i] - 1, nrows)];

      if (debuglevel > 0) {

        if (res[i] == NA_INTEGER)
          Rprintf("  > prediction for observation %d is NA with probabilities:\n");
        else
          Rprintf("  > prediction for observation %d is '%s' with probabilities:\n",
            i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

        Rprintf("  ");
        for (int k = 0; k < nrows; k++)
          Rprintf("  %lf", (cpt + nrows * (configs[i] - 1))[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(1, nmax[configs[i] - 1], res + i, maxima + (configs[i] - 1) * nrows);

      if (debuglevel > 0) {

        Rprintf("  > there are %d levels tied for prediction of observation %d, applying tie breaking.\n",
          nmax[configs[i] - 1], i + 1);
        Rprintf("  > tied levels are:");
        for (k = 0; k < nmax[configs[i] - 1]; k++)
          Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[CMC(k, configs[i] - 1, nrows)] - 1)));
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

    /* attach the prediction probabilities to the return value. */
    if (include_prob) {

      memcpy(pt + i * tr_nlevels, cpt + nrows * (configs[i] - 1),
        tr_nlevels * sizeof(double));

    }/*THEN*/

  }/*FOR*/

  /* save the state of the random number generator. */
  PutRNGstate();

  if (include_prob > 0) {

    /* set the levels of the target variable as rownames. */
    setDimNames(probtab, tr_levels, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(2);

  }/*THEN*/
  else {

    UNPROTECT(1);

  }/*ELSE*/

  Free1D(iscratch);
  Free1D(buf);
  Free1D(dscratch);
  Free1D(maxima);
  Free1D(nmax);

  return result;

}/*CDPRED*/

/* predict the values of a conditional Gaussian node. */
SEXP ccgpred(SEXP fitted, SEXP configurations, SEXP parents, SEXP debug) {

int i = 0, j = 0, ndata = length(configurations);
int *config = INTEGER(configurations), debuglevel = isTRUE(debug);
int cur_config = 0, np = length(parents), nrows = np + 1;
double *res = NULL, *beta = NULL, *beta_offset = NULL, **columns = NULL;
SEXP result;

  /* get the regression coefficients of the conditional Gaussian distribution. */
  beta = REAL(getListElement(fitted, "coefficients"));

  /* dereference the columns of the data frame and get the sample size. */
  columns = (double **) Calloc1D(np, sizeof(double *));
  for (i = 0; i < np; i++)
    columns[i] = REAL(VECTOR_ELT(parents, i));

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(REALSXP, ndata));
  res = REAL(result);

  for (i = 0; i < ndata; i++)  {

    /* find out which conditional regression to use for prediction. */
    cur_config = config[i] - 1;
    beta_offset = beta + cur_config * nrows;

    /* compute the mean value for this observation. */
    res[i] = beta_offset[0];

    for (j = 0; j < np; j++)
      res[i] += columns[j][i] * beta_offset[j + 1];

    if (debuglevel > 0) {

      Rprintf("  > prediction for observation %d is %lf with predictor:\n",
        i + 1, res[i]);

      Rprintf("    (%lf)", beta_offset[0]);
      for (j = 0; j < np; j++)
        Rprintf(" + (%lf) * (%lf)", columns[j][i], beta_offset[j + 1]);
      Rprintf("\n");

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  Free1D(columns);

  return result;

}/*CCGPRED*/

/* predict the value of the training variable in a naive Bayes or Tree-Augmented
 * naive Bayes classifier. */
SEXP naivepred(SEXP fitted, SEXP data, SEXP parents, SEXP training, SEXP prior,
    SEXP prob, SEXP debug) {

int i = 0, j = 0, k = 0, n = 0, nvars = length(fitted), nmax = 0, tr_nlevels = 0;
int *res = NULL, **ex = NULL, *ex_nlevels = NULL;
int idx = 0, *tr_id = INTEGER(training), include_prob = isTRUE(prob);
int *iscratch = NULL, *maxima = NULL, *prn = NULL, debuglevel = isTRUE(debug);
double **cpt = NULL, *pr = NULL, *scratch = NULL, *buf = NULL, *pt = NULL;
double sum = 0;
SEXP temp, tr, tr_levels, tr_class, result, nodes, probtab = R_NilValue;

  /* cache the node labels. */
  nodes = getAttrib(fitted, R_NamesSymbol);

  /* cache the pointers to all the variables. */
  ex = (int **) Calloc1D(nvars, sizeof(int *));
  ex_nlevels = Calloc1D(nvars, sizeof(int));

  for (i = 0; i < nvars; i++) {

    temp = VECTOR_ELT(data, i);
    ex[i] = INTEGER(temp);
    ex_nlevels[i] = NLEVELS(temp);

  }/*FOR*/

  /* get the training variable and its levels. */
  n = length(VECTOR_ELT(data, 0));
  tr = getListElement(VECTOR_ELT(fitted, *tr_id - 1), "prob");
  tr_levels = VECTOR_ELT(getAttrib(tr, R_DimNamesSymbol), 0);
  tr_class = getAttrib(VECTOR_ELT(data, *tr_id - 1), R_ClassSymbol);
  tr_nlevels = length(tr_levels);
  /* get the prior distribution. */
  pr = REAL(prior);

  if (debuglevel > 0) {

    Rprintf("* the prior distribution for the target variable is:\n");
    PrintValue(prior);

  }/*THEN*/

  /* allocate the scratch space used to compute posterior probabilities. */
  scratch = Calloc1D(tr_nlevels, sizeof(double));
  buf = Calloc1D(tr_nlevels, sizeof(double));

  /* cache the pointers to the conditional probability tables. */
  cpt = (double **) Calloc1D(nvars, sizeof(double *));

  for (i = 0; i < nvars; i++)
    cpt[i] = REAL(getListElement(VECTOR_ELT(fitted, i), "prob"));

  /* dereference the parents' vector. */
  prn = INTEGER(parents);

  /* create the vector of indexes. */
  iscratch = Calloc1D(tr_nlevels, sizeof(int));

  /* allocate the array for the indexes of the maxima. */
  maxima = Calloc1D(tr_nlevels, sizeof(int));

  /* allocate the return value. */
  PROTECT(result = allocVector(INTSXP, n));
  res = INTEGER(result);

  /* allocate and initialize the table of the posterior probabilities. */
  if (include_prob > 0) {

    PROTECT(probtab = allocMatrix(REALSXP, tr_nlevels, n));
    pt = REAL(probtab);
    memset(pt, '\0', n * tr_nlevels * sizeof(double));

  }/*THEN*/

  /* initialize the random seed, just in case we need it for tie breaking. */
  GetRNGstate();

  /* for each observation... */
  for (i = 0; i < n; i++) {

    /* ... reset the scratch space and the indexes array... */
    for (k = 0; k < tr_nlevels; k++) {

      scratch[k] = log(pr[k]);
      iscratch[k] = k + 1;

    }/*FOR*/

    if (debuglevel > 0)
      Rprintf("* predicting the value of observation %d.\n", i + 1);

    /* ... and for each conditional probability table... */
    for (j = 0; j < nvars; j++) {

      /* ... skip the training variable... */
      if (*tr_id == j + 1)
        continue;

      /* ... (this is the root node of the Chow-Liu tree) ... */
      if (prn[j] == NA_INTEGER) {

        /* ... and for each row of the conditional probability table... */
        for (k = 0; k < tr_nlevels; k++) {

          if (debuglevel > 0) {

            Rprintf("  > node %s: picking cell %d (%d, %d) from the CPT (p = %lf).\n",
              NODE(j), CMC(ex[j][i] - 1, k, ex_nlevels[j]), ex[j][i], k + 1,
              cpt[j][CMC(ex[j][i] - 1, k, ex_nlevels[j])]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] += log(cpt[j][CMC(ex[j][i] - 1, k, ex_nlevels[j])]);

        }/*FOR*/

      }/*THEN*/
      else {

        /* ... and for each row of the conditional probability table... */
        for (k = 0; k < tr_nlevels; k++) {

          /* (the first dimension corresponds to the current node [X], the second
           * to the training node [Y], the third to the only parent of the current
           * node [Z]; CMC coordinates are computed as X + Y * NX + Z * NX * NY. */
          idx = (ex[j][i] - 1) + k * ex_nlevels[j] +
                  (ex[prn[j] - 1][i] - 1) * ex_nlevels[j] * tr_nlevels;

          if (debuglevel > 0) {

            Rprintf("  > node %s: picking cell %d (%d, %d, %d) from the CPT (p = %lf).\n",
              NODE(j), idx, ex[j][i], k + 1, ex[prn[j] - 1][i], cpt[j][idx]);

          }/*THEN*/

          /* ... update the posterior probability. */
          scratch[k] += log(cpt[j][idx]);

        }/*FOR*/

      }/*ELSE*/

    }/*FOR*/

    /* find out the mode(s). */
    nmax = all_max(scratch, tr_nlevels, maxima, iscratch, buf);

    /* compute the posterior probabilities on the right scale, to attach them
     * to the return value. */
    if (include_prob) {

      /* copy the log-probabilities from scratch. */
      memcpy(pt + i * tr_nlevels, scratch, tr_nlevels * sizeof(double));

      /* transform log-probabilities into plain probabilities. */
      for (k = 0, sum = 0; k < tr_nlevels; k++)
        sum += pt[i * tr_nlevels + k] = exp(pt[i * tr_nlevels + k] - scratch[maxima[0] - 1]);

      /* rescale them to sum up to 1. */
      for (k = 0; k < tr_nlevels; k++)
        pt[i * tr_nlevels + k] /= sum;

    }/*THEN*/

    if (nmax == 1) {

      res[i] = maxima[0];

      if (debuglevel > 0) {

        Rprintf("  @ prediction for observation %d is '%s' with (log-)posterior:\n",
          i + 1, CHAR(STRING_ELT(tr_levels, res[i] - 1)));

        Rprintf("  ");
        for (k = 0; k < tr_nlevels; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

      }/*THEN*/

    }/*THEN*/
    else {

      /* break ties: sample with replacement from all the maxima. */
      SampleReplace(1, nmax, res + i, maxima);

      if (debuglevel > 0) {

        Rprintf("  @ there are %d levels tied for prediction of observation %d, applying tie breaking.\n", nmax, i + 1);

        Rprintf("  ");
        for (k = 0; k < tr_nlevels; k++)
          Rprintf("  %lf", scratch[k]);
        Rprintf("\n");

        Rprintf("  @ tied levels are:");
        for (k = 0; k < nmax; k++)
          Rprintf(" %s", CHAR(STRING_ELT(tr_levels, maxima[k] - 1)));
        Rprintf(".\n");

      }/*THEN*/

    }/*ELSE*/

  }/*FOR*/

  /* save the state of the random number generator. */
  PutRNGstate();

  /* add back the attributes and the class to the return value. */
  setAttrib(result, R_LevelsSymbol, tr_levels);
  setAttrib(result, R_ClassSymbol, tr_class);

  if (include_prob > 0) {

    /* set the levels of the target variable as rownames. */
    setDimNames(probtab, tr_levels, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(2);

  }/*THEN*/
  else {

    UNPROTECT(1);

  }/*ELSE*/

  Free1D(ex);
  Free1D(ex_nlevels);
  Free1D(scratch);
  Free1D(buf);
  Free1D(cpt);
  Free1D(iscratch);
  Free1D(maxima);

  return result;

}/*NAIVEPRED*/

