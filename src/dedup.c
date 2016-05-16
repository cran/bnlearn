#include "include/rcore.h"
#include "include/dataframe.h"
#include "include/globals.h"
#include "include/covariance.h"

/* remove one variable in each highly-correlated pair. */
SEXP dedup (SEXP data, SEXP threshold, SEXP debug) {

int i = 0, j = 0, k = 0, dropped = 0;
int ncols = length(data), nrows = length(VECTOR_ELT(data, 0));
int debuglevel = isTRUE(debug);
short int *drop = NULL;
double **column = NULL, *mean = NULL, *sse = NULL;
double tol = MACHINE_TOL, t = NUM(threshold);
long double sum = 0;
SEXP result, colnames, nodes = getAttrib(data, R_NamesSymbol);

  /* set up a counter to flag the variables. */
  drop = Calloc1D(ncols, sizeof(short int));

  /* extract the columns from the data frame. */
  column = Calloc1D(ncols, sizeof(double *));
  for (j = 0; j < ncols; j++)
    column[j] = REAL(VECTOR_ELT(data, j));

  if (debuglevel > 0)
    Rprintf("* caching means.\n");

  /* cache the means. */
  mean = Calloc1D(ncols, sizeof(double));
  c_meanvec(column, mean, nrows, ncols, 0);

  if (debuglevel > 0)
    Rprintf("* caching standard deviations.\n");

  /* cache the variances. */
  sse = Calloc1D(ncols, sizeof(double));
  c_ssevec(column, sse, mean, nrows, ncols, 0);

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

      /* safety check against "divide by zero" errors. */
      if ((sse[j] < tol) || (sse[k] < tol))
        sum = 0;
      else
        sum /= sqrt(sse[j] * sse[k]);

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

  Free1D(drop);
  Free1D(column);
  Free1D(mean);
  Free1D(sse);

  UNPROTECT(2);

  return result;

}/*DEDUP*/

