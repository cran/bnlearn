
#include <Rmath.h>
#include "common.h"

/* posterior dirichlet probability (covers BDe and K2 scores). */
SEXP dpost (SEXP x, SEXP lx, SEXP length, SEXP iss, SEXP debug) {

int k = 0;
int *llx = INTEGER(lx), *num = INTEGER(length), *xx = INTEGER(x);
int *n = NULL, *imaginary = NULL;
double alpha = 0, *res = NULL;
SEXP result;

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distriburion. */
    imaginary = llx;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INTEGER(iss);
    alpha = (double) *imaginary / (double) *llx;

  }/*ELSE*/

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table. */
  n = alloc1dcont(*llx);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++)
    n[xx[k] - 1]++;

  if (isTRUE(debug)) {

    Rprintf("  > prior distribution is: %lf (x%d)\n", alpha, *llx);
    Rprintf("  > sample distribution is:\n    ");
    for (k = 0; k < *llx; k++) {

      Rprintf("%d ", n[k]);
      if (((k % 6) == 0) && (k != 0)) Rprintf("\n    ");

    }/*FOR*/
    Rprintf("\n");
    Rprintf("  > posterior distribution is:\n    ");
    for (k = 0; k < *llx; k++) {

      Rprintf("%lf ", n[k] + alpha);
      if (((k % 6) == 0) && (k != 0)) Rprintf("\n    ");

    }/*FOR*/
    Rprintf("\n");
    Rprintf("  > real sample size: %d\n", *num);

  }/*THEN*/

  /* compute the posterior probability. */
  for (k = 0; k < *llx; k++)
    *res += lgammafn(n[k] + alpha) - lgammafn(alpha);
  *res += lgammafn((double)(*imaginary)) -
            lgammafn((double)(*imaginary + *num));

  UNPROTECT(1);
  return result;

}/*DPOST*/

/* conditional posterior dirichlet probability (covers BDe and K2 scores). */
SEXP cdpost (SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length, SEXP iss,
    SEXP nparams, SEXP debug) {

int j = 0, k = 0;
int *llx = INTEGER(lx), *lly = INTEGER(ly), *p = INTEGER(nparams);
int *xx = INTEGER(x), *yy = INTEGER(y), *num = INTEGER(length);
int *imaginary = NULL, **n = NULL, *nj = NULL;
double alpha = 0, *res = NULL;
SEXP result;

  if (isNull(iss)) {

    /* this is for K2, which does not define an imaginary sample size;
     * all hyperparameters are set to 1 in the prior distriburion. */
    imaginary = p;
    alpha = 1;

  }/*THEN*/
  else {

    /* this is for the BDe score. */
    imaginary = INTEGER(iss);
    alpha = (double) *imaginary / (double) *p;

  }/*ELSE*/

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  res = REAL(result);
  *res = 0;

  /* initialize the contingency table. */
  n = alloc2dcont(*llx, *lly);
  nj = alloc1dcont(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;
    nj[yy[k] - 1]++;

  }/*FOR*/

  if (isTRUE(debug)) {

    Rprintf("  > prior distribution is: %lf (x%d)\n", alpha, *llx);
    Rprintf("  > sample distribution is:\n");
    for (j = 0; j < *llx; j++) {

      Rprintf("    ");

      for (k = 0; k < *lly; k++) {

        Rprintf("%d ", n[j][k]);

      }/*FOR*/

      Rprintf("\n");

    }/*FOR*/
    Rprintf("  > posterior distribution is:\n");
    for (j = 0; j < *llx; j++) {

      Rprintf("    ");

      for (k = 0; k < *lly; k++) {

        Rprintf("%lf ", n[j][k] + alpha);

      }/*FOR*/

      Rprintf("\n");

    }/*FOR*/

    Rprintf("  > real sample size: %d\n", *num);

  }/*THEN*/

  /* compute the conditional posterior probability. */
  for (j = 0; j < *llx; j++)
    for (k = 0; k < *lly; k++)
      *res += lgammafn(n[j][k] + alpha) - lgammafn(alpha);
  for (k = 0; k < *lly; k++)
    *res += lgammafn((double)*imaginary / *lly) -
              lgammafn(nj[k] + (double)*imaginary / *lly);

  UNPROTECT(1);
  return result;

}/*CDPOST*/
