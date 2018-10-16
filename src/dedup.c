#include "include/rcore.h"
#include "include/data.frame.h"
#include "include/globals.h"
#include "include/covariance.h"
#include "include/data.table.h"

/* remove one variable in each highly-correlated pair. */
SEXP dedup (SEXP data, SEXP threshold, SEXP complete, SEXP debug) {

int i = 0, j = 0, k = 0, dropped = 0, nc = 0;
int debuglevel = isTRUE(debug);
double *mean = NULL, *sse = NULL, *xx = NULL, *yy = NULL;
double cur_mean[2], cur_sse[2];
double tol = MACHINE_TOL, t = NUM(threshold);
long double sum = 0;
SEXP result, colnames;
gdata dt = { 0 };

  /* extract the columns from the data frame. */
  dt = gdata_from_SEXP(data, 0);
  meta_init_flags(&(dt.m), 0, complete, R_NilValue);
  meta_copy_names(&(dt.m), 0, data);
  /* set up the vectors for the pairwise complete observations. */
  xx = Calloc1D(dt.m.nobs, sizeof(double));
  yy = Calloc1D(dt.m.nobs, sizeof(double));

  if (debuglevel > 0)
    Rprintf("* caching means and variances.\n");

  mean = Calloc1D(dt.m.ncols, sizeof(double));
  sse = Calloc1D(dt.m.ncols, sizeof(double));

  /* cache the mean and variance of complete variables. */
  for (j = 0; j < dt.m.ncols; j++) {

    if (!dt.m.flag[j].complete)
      continue;

    mean[j] = c_mean(dt.col[j], dt.m.nobs);
    sse[j] = c_sse(dt.col[j], mean[j], dt.m.nobs);

  }/*FOR*/

  /* main loop. */
  for (j = 0; j < dt.m.ncols - 1; j++) {

    /* skip variables already flagged for removal. */
    if (dt.m.flag[j].drop)
      continue;

    if (debuglevel > 0)
      Rprintf("* looking at %s with %d variables still to check.\n",
        dt.m.names[j], dt.m.ncols - (j + 1));

    for (k = j + 1; k < dt.m.ncols; k++) {

      /* skip variables already flagged for removal. */
      if (dt.m.flag[k].drop)
        continue;

      if (dt.m.flag[j].complete && dt.m.flag[k].complete) {

        /* use the cached means and variances. */
        cur_mean[0] = mean[j];
        cur_mean[1] = mean[k];
        cur_sse[0] = sse[j];
        cur_sse[1] = sse[k];

        /* compute the covariance. */
        for (i = 0, sum = 0; i < dt.m.nobs; i++)
          sum += (dt.col[j][i] - cur_mean[0]) * (dt.col[k][i] - cur_mean[1]);

      }/*THEN*/
      else {

        for (i = 0, nc = 0; i < dt.m.nobs; i++) {

          if (ISNAN(dt.col[j][i]) || ISNAN(dt.col[k][i]))
            continue;

          xx[nc] = dt.col[j][i];
          yy[nc++] = dt.col[k][i];

        }/*FOR*/


        /* if there are no complete observations, take the variables to be
         * independent. */
        if (nc == 0)
          continue;

        cur_mean[0] = c_mean(xx, nc);
        cur_mean[1] = c_mean(yy, nc);
        cur_sse[0] = c_sse(xx, cur_mean[0], nc);
        cur_sse[1] = c_sse(yy, cur_mean[1], nc);

        /* compute the covariance. */
        for (i = 0, sum = 0; i < nc; i++)
          sum += (xx[i] - cur_mean[0]) * (yy[i] - cur_mean[1]);

      }/*ELSE*/

      /* safety check against "divide by zero" errors. */
      if ((cur_sse[0] < tol) || (cur_sse[1] < tol))
        sum = 0;
      else
        sum /= sqrt(cur_sse[0] * cur_sse[1]);

      /* test the correlation against the threshold. */
      if (fabsl(sum) > t) {

        if (debuglevel > 0)
          Rprintf("%s is collinear with %s, dropping %s with COR = %.4Lf\n",
            dt.m.names[j], dt.m.names[k], dt.m.names[k], sum);

        /* flag the variable for removal. */
        dt.m.flag[k].drop = TRUE;
        dropped++;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  /* set up the return value. */
  PROTECT(result = allocVector(VECSXP, dt.m.ncols - dropped));
  PROTECT(colnames = allocVector(STRSXP, dt.m.ncols - dropped));

  for (j = 0, k = 0; j < dt.m.ncols; j++)
    if (!dt.m.flag[j].drop) {

      SET_STRING_ELT(colnames, k, mkChar(dt.m.names[j]));
      SET_VECTOR_ELT(result, k++, VECTOR_ELT(data, j));

    }/*THEN*/

  setAttrib(result, R_NamesSymbol, colnames);

  /* make it a data frame. */
  minimal_data_frame(result);

  Free1D(mean);
  Free1D(sse);
  Free1D(xx);
  Free1D(yy);
  FreeGDT(dt, FALSE);

  UNPROTECT(2);

  return result;

}/*DEDUP*/

