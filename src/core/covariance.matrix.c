#include "../include/rcore.h"
#include "allocations.h"
#include "../math/linear.algebra.h"
#include "covariance.matrix.h"

covariance new_covariance(int dim, bool decomp) {

covariance cov = { 0 };

  cov.mat = Calloc1D(dim * dim, sizeof(double));
  cov.dim = dim;

  if (decomp) {

    cov.u = Calloc1D(dim * dim, sizeof(double));
    cov.d = Calloc1D(dim, sizeof(double));
    cov.vt = Calloc1D(dim * dim, sizeof(double));

  }/*THEN*/

  return cov;

}/*NEW_COVARIANCE*/

void print_covariance(covariance cov) {

  for (int i = 0; i < cov.dim; i++) {

    for (int j = 0; j < cov.dim; j++)
      Rprintf("%lf ", cov.mat[CMC(i, j, cov.dim)]);

   Rprintf("\n");

  }/*FOR*/

}/*PRINT_COVARIANCE*/

void copy_covariance(covariance *src, covariance *copy) {

  memcpy((*copy).mat, (*src).mat, (*copy).dim * (*copy).dim * sizeof(double));

}/*COPY_COVARIANCE*/

void covariance_drop_variable(covariance *full, covariance *sub, int to_drop) {

  for (int j = 0, k = 0; j < (*full).dim; j++)
    for (int i = 0; i < (*full).dim; i++)
      if ((i != to_drop) && (j != to_drop))
        (*sub).mat[k++] = (*full).mat[CMC(i, j, (*full).dim)];

  (*sub).dim = (*full).dim - 1;

}/*COVARIANCE_DROP_COLUMN*/

void FreeCOV(covariance cov) {

  Free1D(cov.mat);
  Free1D(cov.u);
  Free1D(cov.d);
  Free1D(cov.vt);

}/*FREECOV*/

/* fill a covariance matrix. */
void c_covmat(double **data, double *mean, int nrow, int ncol,
    covariance cov, int first) {

int i = 0, j = 0, k = 0;
long double temp = 0;

  /* special case: with zero and one observations, the covariance
   * matrix is all zeroes. */
  if (nrow <= 1) {

    memset(cov.mat, '\0', ncol * ncol * sizeof(double));

    return ;

  }/*THEN*/

  /* compute the actual covariance. */
  for (i = first; i < ncol; i++) {

    for (j = i; j < ncol; j++) {

      for (k = 0, temp = 0; k < nrow; k++)
        temp += (data[i][k] - mean[i]) * (data[j][k] - mean[j]);

      /* fill in the symmetric element of the matrix. */
      cov.mat[CMC(j, i, ncol)] = cov.mat[CMC(i, j, ncol)] =
        (double)(temp / (nrow - 1));

    }/*FOR*/

  }/*FOR*/

}/*C_COVMAT*/

/* fill a covariance matrix from incomplete data, computing means as well. */
void c_covmat_with_missing(double **data, int nrow, int ncol,
    bool *missing_partial, bool *missing_all, double *mean, double *mat,
    int *ncomplete) {

int i = 0, j = 0, k = 0, nc = 0;
long double temp = 0;

  /* clear the overall missingness indicators, since this function is always
   * called in loops with allocated-once scratch spaces. */
  memset(missing_all, '\0', nrow * sizeof(bool));

  /* count the number of complete observations. */
  for (i = 0; i < nrow; i++) {

    if (missing_partial) {

      if (missing_partial[i]) {

        missing_all[i] = TRUE;
        continue;

      }/*THEN*/

    }/*THEN*/

    for (j = 0; j < ncol; j++) {

      if (ISNAN(data[j][i])) {

        missing_all[i] = TRUE;
        break;

      }/*THEN*/

    }/*FOR*/

    if (!missing_all[i])
      nc++;

  }/*FOR*/

  *ncomplete = nc;

  /* if there are no complete data points, return a zero matrix. */
  if (nc == 0)
    return;

  /* compute the means. */
  for (j = 0; j < ncol; j++) {

    for (i = 0, temp = 0; i < nrow; i++) {

      if (missing_all[i])
        continue;

      temp += data[j][i];

    }/*FOR*/

    mean[j] = temp / nc;

  }/*FOR*/

  /* compute the covariance from complete observations. */
  for (i = 0; i < ncol; i++) {
    for (j = 0; j < ncol; j++) {

      for (k = 0, temp = 0; k < nrow; k++) {

        if (missing_all[k])
          continue;

        temp += (data[j][k] - mean[j]) * (data[i][k] - mean[i]);

      }/*FOR*/

      /* fill in the symmetric element of the matrix. */
      mat[CMC(j, i, ncol)] = mat[CMC(i, j, ncol)] =
        (double)(temp / (nc - 1));

    }/*FOR*/

  }/*FOR*/

}/*C_COVMAT_WITH_MISSING*/

/* update only a single row/column in a covariance matrix. */
void c_update_covmat(double **data, double *mean, int update, int nrow,
    int ncol, double *mat) {

int j = 0, k = 0;
long double temp = 0;

  /* special case: with zero and one observations, set the row/column to zero. */
  if (nrow <= 1) {

    for (j = 0; j < ncol; j++)
      mat[CMC(j, update, ncol)] = mat[CMC(update, j, ncol)] = 0;

    return ;

  }/*THEN*/

  /* compute the actual covariance. */
  for (j = 0; j < ncol; j++) {

    for (k = 0, temp = 0; k < nrow; k++)
      temp += (data[update][k] - mean[update]) * (data[j][k] - mean[j]);

    /* fill the symmetric elements of the matrix. */
    mat[CMC(j, update, ncol)] = mat[CMC(update, j, ncol)] =
      (double)(temp / (nrow - 1));

  }/*FOR*/

}/*C_UPDATE_COVMAT*/

