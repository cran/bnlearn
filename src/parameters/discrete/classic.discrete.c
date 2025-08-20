#include "../../include/rcore.h"
#include "../../core/contingency.tables.h"
#include "../../math/linear.algebra.h"

void c_classic_discrete_parameters(int *counts, double *cpt, int nrows,
    int ncols, double alpha, bool replace) {

long double colsum = 0;

  /* add the imaginary sample size, if any, ... */
  for (int i = 0; i < nrows * ncols; i++)
    cpt[i] = counts[i] + alpha / (nrows * ncols);

  /* ... and normalize the columns to sum up to 1. */
  for (int j = 0; j < ncols; j++) {

    colsum = 0;
    for (int i = 0; i < nrows; i++)
      colsum += cpt[CMC(i, j, nrows)];

    /* some columns cannot be normalized: either fill them with NaNs or with
     * a uniform distribution. */
    if (colsum == 0) {

      if (replace) {

        for (int i = 0; i < nrows; i++)
          cpt[CMC(i, j, nrows)] = (double) 1 / nrows;

      }/*THEN*/
      else {

        for (int i = 0; i < nrows; i++)
          cpt[CMC(i, j, nrows)] = NA_REAL;

      }/*ELSE*/

    }/*THEN*/
    else {

      for (int i = 0; i < nrows; i++)
        cpt[CMC(i, j, nrows)] /= colsum;

    }/*ELSE*/

  }/*FOR*/

}/*C_CLASSIC_DISCRETE_PARAMETERS*/

