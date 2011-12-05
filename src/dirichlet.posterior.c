
#include <Rmath.h>
#include "common.h"

/* posterior dirichlet probability (covers BDe and K2 scores). */
SEXP dpost(SEXP x, SEXP iss, SEXP exp) {

int i = 0, k = 0, num = LENGTH(x);
int llx = NLEVELS(x), *xx = INTEGER(x), *n = NULL, *imaginary = NULL;
double alpha = 0, *res = NULL;
SEXP result;

  /* the correct vaules for the hyperparameters alpha are documented in
   * "Learning Bayesian Networks: The Combination of Knowledge and Statistical
   * Data" by Heckerman, Geiger & Chickering (1995), page 17. */

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distribution. */
    imaginary = &llx;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INTEGER(iss);
    alpha = (double) *imaginary / (double) llx;

  }/*ELSE*/

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table. */
  n = alloc1dcont(llx);

  /* compute the frequency table of x, disregarding experimental data. */
  if (exp == R_NilValue) {

    for (i = 0; i < num; i++)
      n[xx[i] - 1]++;

  }/*THEN*/
  else {

    int *e = INTEGER(exp);

    for (i = 0, k = 0; i < num; i++) {
      if (i != e[k] - 1)
        n[xx[i] - 1]++;
      else
        k++;

    }/*FOR*/

    /* adjust the sample size to match the number of observational data. */
   num -= LENGTH(exp);

  }/*ELSE*/

  /* compute the posterior probability. */
  for (i = 0; i < llx; i++)
    *res += lgammafn(n[i] + alpha) - lgammafn(alpha);
  *res += lgammafn((double)(*imaginary)) -
            lgammafn((double)(*imaginary + num));

  UNPROTECT(1);
  return result;

}/*DPOST*/

/* conditional posterior dirichlet probability (covers BDe and K2 scores). */
SEXP cdpost(SEXP x, SEXP y, SEXP iss, SEXP exp, SEXP nparams) {

int i = 0, j = 0, k = 0, imaginary = 0, num = LENGTH(x);
int llx = NLEVELS(x), lly = NLEVELS(y), *xx = INTEGER(x), *yy = INTEGER(y);
int **n = NULL, *nj = NULL;
double alpha = 0, *res = NULL, *p = REAL(nparams);
SEXP result;

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distribution. */
    imaginary = (int) *p;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INT(iss);
    alpha = (double) imaginary / *p;

  }/*ELSE*/

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table. */
  n = alloc2dcont(llx, lly);
  nj = alloc1dcont(lly);

  /* compute the joint frequency of x and y. */
  if (exp == R_NilValue) {

    for (i = 0; i < num; i++) {

      n[xx[i] - 1][yy[i] - 1]++;
      nj[yy[i] - 1]++;

    }/*FOR*/

  }/*THEN*/
  else {

    int *e = INTEGER(exp);

    for (i = 0, k = 0; i < num; i++) {

      if (i != e[k] - 1) {

        n[xx[i] - 1][yy[i] - 1]++;
        nj[yy[i] - 1]++;

      }/*THEN*/
      else {

        k++;

      }/*ELSE*/

    }/*FOR*/

    /* adjust the sample size to match the number of observational data. */
   num -= LENGTH(exp);

  }/*ELSE*/

  /* compute the conditional posterior probability. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      *res += lgammafn(n[i][j] + alpha) - lgammafn(alpha);
  for (j = 0; j < lly; j++)
    *res += lgammafn((double)imaginary / lly) -
              lgammafn(nj[j] + (double)imaginary / lly);

  UNPROTECT(1);
  return result;

}/*CDPOST*/
