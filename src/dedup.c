#include "include/rcore.h"
#include "include/dataframe.h"
#include "include/allocations.h"
#include "include/globals.h"

SEXP dedup (SEXP data, SEXP threshold, SEXP debug) {

int i = 0, j = 0, k = 0, dropped = 0;
int ncols = length(data), nrows = length(VECTOR_ELT(data, 0));
int debuglevel = isTRUE(debug);
short int *drop = NULL;
double **column = NULL, *mean = NULL, *sd = NULL;
double tol = MACHINE_TOL, t = NUM(threshold);
long double sum = 0;
SEXP result, colnames, nodes = getAttrib(data, R_NamesSymbol);

  /* set up a counter to flag the variables. */
  drop = allocstatus(ncols);

  /* extract the columns from the data frame. */
  column = (double **) alloc1dpointer(ncols);
  for (j = 0; j < ncols; j++)
    column[j] = REAL(VECTOR_ELT(data, j));

  if (debuglevel > 0)
    Rprintf("* caching means.\n");

  /* cache the means. */
  mean = alloc1dreal(ncols);
  for (j = 0; j < ncols; j++) {

    for (i = 0; i < nrows; i++)
      mean[j] += column[j][i];
    mean[j] /= nrows;

  }/*FOR*/

  if (debuglevel > 0)
    Rprintf("* caching standard deviations.\n");

  /* cache the variances. */
  sd = alloc1dreal(ncols);
  for (j = 0; j < ncols; j++) {

    for (i = 0, sum = 0; i < nrows; i++)
      sum += (column[j][i] - mean[j]) * (column[j][i] - mean[j]);
    sd[j] = sqrt(sum / nrows);

  }/*FOR*/

  /* main loop. */
  for (j = 0; j < ncols - 1; j++) {

    /* skip variables already flagged for removal. */
    if (drop[j])
      continue;

    if (debuglevel > 0)
      Rprintf("* looking at %s with %d variables still to check.\n",
        NODE(j), ncols - (j + 1));

    for (k = j + 1; k < ncols; k++) {

      /* skip variables already flagged for removal. */
      if (drop[k])
        continue;

      /* compute the covariance. */
      for (i = 0, sum = 0; i < nrows; i++)
        sum += (column[j][i] - mean[j]) * (column[k][i] - mean[k]);
      sum /= nrows;

      /* safety check against "divide by zero" errors. */
      if ((sd[j] < tol) || (sd[k] < tol))
        sum = 0;
      else
        sum /= sd[j] * sd[k];

      /* test the correlation against the threshold. */
      if (fabsl(sum) > t) {

        if (debuglevel > 0)
          Rprintf("%s is collinear with %s, dropping %s with COR = %.4Lf\n",
            NODE(j), NODE(k), NODE(k), sum);

        /* flag the variable for removal. */
        drop[k] = TRUE;
        dropped++;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  /* set up the return value. */
  PROTECT(result = allocVector(VECSXP, ncols - dropped));
  PROTECT(colnames = allocVector(STRSXP, ncols - dropped));

  for (j = 0, k = 0; j < ncols; j++)
    if (!drop[j]) {

      SET_STRING_ELT(colnames, k, STRING_ELT(nodes, j));
      SET_VECTOR_ELT(result, k++, VECTOR_ELT(data, j));

    }/*THEN*/

  setAttrib(result, R_NamesSymbol, colnames);

  /* make it a data frame. */
  minimal_data_frame(result);

  UNPROTECT(2);

  return result;

}/*DEDUP*/

