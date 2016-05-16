#include "include/rcore.h"
#include "include/dataframe.h"
#include "include/sampling.h"
#include "include/globals.h"
#include "include/matrix.h"

static double posterior_mean(double *x, double *wgt, int n, int debuglevel) {

int k = 0;
long double wsum = 0, wtot = 0;

  for (k = 0; k < n; k++) {

    wsum += x[k] * wgt[k];
    wtot += wgt[k];

  }/*FOR*/
  wsum /= wtot;

  /* if all weights are zero, the predicted value is NA.*/
  if (wtot == 0)
    wsum = NA_REAL;

  if (debuglevel > 0)
    Rprintf("  > prediction is %Lf.\n", wsum);

  return (double)wsum;

}/*POSTERIOR_MEAN*/

static int posterior_mode(int *x, double *wgt, int n, long double *counts,
    SEXP levels, int nlvls, int *drop, int debuglevel) {

int k = 0, res = 0;

  memset(counts, '\0', nlvls * sizeof(long double));

  for (k = 0; k < n; k++) {

    /* c_rbn_master() may generate NAs, disregard and print a warning. */
    if (x[k] == NA_INTEGER)
      (*drop)++;
    else
      counts[x[k] - 1] += wgt[k];

  }/*FOR*/

  res = ld_which_max(counts, nlvls);

  /* if all weights are zero, the predicted value is NA.*/
  if (counts[res - 1] == 0)
    res = NA_INTEGER;

  if (debuglevel > 0) {

    Rprintf("  > prediction is '%s' with weight sums:\n",
      res == NA_INTEGER ? "NA" : CHAR(STRING_ELT(levels, res - 1)));
    for (k = 0; k < nlvls; k++)
      Rprintf("%Lf ", counts[k]);
    Rprintf("\n");

  }/*THEN*/

  return res;

}/*POSTERIOR_MODE*/

/* predict the values of one or more variables given one or more variables by
 * maximum a posteriori (MAP). */
SEXP mappred(SEXP node, SEXP fitted, SEXP data, SEXP n, SEXP from, SEXP prob,
    SEXP debug) {

int i = 0, j = 0, nobs = 0, nev = 0, nlvls = 0, drop = 0;
int *vartypes = NULL, nsims = INT(n), debuglevel = isTRUE(debug);
int include_prob = INT(prob);
void **varptrs = NULL, **evptrs = NULL, *pred = NULL, *res = NULL;
SEXP result, colnames, evidence, evmatch, temp = R_NilValue;
SEXP cpdist, predicted, lvls = R_NilValue, probtab = R_NilValue;
double *wgt = NULL, *pt = NULL;
long double *lvls_counts = NULL, lvls_tot = 0;

  /* extract the names of the variables in the data. */
  colnames = getAttrib(data, R_NamesSymbol);

  /* remove the name of the variable to predict. */
  nev = length(from);
  PROTECT(evmatch = match(colnames, from, 0));

  /* cache variable types and pointers. */
  vartypes = Calloc1D(nev, sizeof(int));
  varptrs = (void **) Calloc1D(nev, sizeof(void *));
  for (j = 0; j < nev; j++) {

    temp = VECTOR_ELT(data, INTEGER(evmatch)[j] - 1);
    vartypes[j] = TYPEOF(temp);
    varptrs[j] = DATAPTR(temp);

  }/*FOR*/

  /* cache the sample size. */
  nobs = length(temp);

  /* allocate a list to hold the evidence. */
  PROTECT(evidence = allocVector(VECSXP, nev));
  setAttrib(evidence, R_NamesSymbol, from);

  /* cache pointers to the elements of the evidence .*/
  evptrs = (void **) Calloc1D(nev, sizeof(void *));

  for (j = 0; j < nev; j++) {

    PROTECT(temp = allocVector(vartypes[j], 1));
    evptrs[j] = DATAPTR(temp);
    SET_VECTOR_ELT(evidence, j, temp);
    UNPROTECT(1);

  }/*FOR*/

  /* make the evidence a data frame to compact debugging output. */
  minimal_data_frame(evidence);

  /* allocate the return value. */
  PROTECT(result = fitnode2df(fitted, STRING_ELT(node, 0), nobs));
  res = DATAPTR(result);

  /* in the case of discrete variables, allocate scratch space for levels'
   * frequencies. */
  if (TYPEOF(result) == INTSXP) {

    lvls = getAttrib(result, R_LevelsSymbol);
    nlvls = length(lvls);
    lvls_counts = Calloc1D(nlvls, sizeof(long double));

    if (include_prob > 0) {

      PROTECT(probtab = allocMatrix(REALSXP, nlvls, nobs));
      pt = REAL(probtab);
      memset(pt, '\0', nobs * nlvls * sizeof(double));

    }/*THEN*/

  }/*THEN*/

  /* allocate the weights. */
  wgt = Calloc1D(nsims, sizeof(double));

  /* allocate sratch space for the random samplings. */
  PROTECT(cpdist = fit2df(fitted, nsims));
  predicted = getListElement(cpdist, (char *)CHAR(STRING_ELT(node, 0)));
  pred = DATAPTR(predicted);

  /* iterate over the observations. */
  for (i = 0; i < nobs; i++) {

    /* copy the values into the list. */
    for (j = 0; j < nev; j++) {

      switch(vartypes[j]) {

        case REALSXP:

          *((double *)evptrs[j]) = ((double *)varptrs[j])[i];
          break;

        case INTSXP:

          *((int *)evptrs[j]) = ((int *)varptrs[j])[i];
          break;

      }/*SWITCH*/

    }/*FOR*/

    if (debuglevel > 0) {

      Rprintf("* predicting observation %d conditional on:\n", i);
      PrintValue(evidence);

    }/*THEN*/

    /* generate samples from the conditional posterior distribution. */
    c_rbn_master(fitted, cpdist, n, evidence, FALSE);
    /* compute the weights. */
    c_lw_weights(fitted, cpdist, nsims, wgt, from, FALSE);

    /* compute the posterior estimate. */
    switch(TYPEOF(predicted)) {

      case REALSXP:

        /* average the predicted values. */
        ((double *)res)[i] = posterior_mean((double *)pred, wgt, nsims,
                               debuglevel);
        break;

      case INTSXP:

        /* pick the most frequent value. */
        ((int *)res)[i] = posterior_mode((int *)pred, wgt, nsims, lvls_counts,
                            lvls, nlvls, &drop, debuglevel);

        /* compute the posterior probabilities on the right scale, to attach
         * them to the return value. */
        if (include_prob > 0) {

          for (j = 0, lvls_tot = 0; j < nlvls; j++) {

            pt[CMC(j, i, nlvls)] = lvls_counts[j];
            lvls_tot += lvls_counts[j];

          }/*FOR*/

          for (j = 0; j < nlvls; j++)
            pt[CMC(j, i, nlvls)] /= lvls_tot;

        }/*THEN*/

        break;

    }/*SWITCH*/

  }/*FOR*/

  if (include_prob > 0) {

    /* set the levels of the taregt variable as rownames. */
    setDimNames(probtab, lvls, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(5);

  }/*THEN*/
  else {

    UNPROTECT(4);

  }/*ELSE*/

  Free1D(vartypes);
  Free1D(varptrs);
  Free1D(evptrs);
  Free1D(wgt);
  if (TYPEOF(result) == INTSXP)
    Free1D(lvls_counts);

  if (drop > 0)
    warning("dropping %d observations because generated samples are NAs.", drop);

  return result;

}/*MAPPRED*/

