#include "include/rcore.h"
#include "include/matrix.h"

/* fill a covariance matrix. */
void c_covmat(double **data, double *mean, int ncols, int nrows,
    double *mat, int first) {

int i = 0, j = 0, k = 0;
long double temp = 0;

  /* compute the actual covariance. */
  for (i = first; i < ncols; i++) {

    for (j = i; j < ncols; j++) {

      for (k = 0, temp = 0; k < nrows; k++)
        temp += (data[i][k] - mean[i]) * (data[j][k] - mean[j]);

      /* fill in the symmetric element of the matrix. */
      mat[CMC(j, i, ncols)] = mat[CMC(i, j, ncols)] =
        (double)(temp / (nrows - 1));

    }/*FOR*/

  }/*FOR*/

}/*C_COVMAT*/

/* update only a single row/column in a covariance matrix. */
void c_update_covmat(double **data, double *mean, int update, int ncols,
    int nrows, double *mat) {

int j = 0, k = 0;
long double temp = 0;

  /* compute the actual covariance. */
  for (j = 0; j < ncols; j++) {

    for (k = 0, temp = 0; k < nrows; k++)
      temp += (data[update][k] - mean[update]) * (data[j][k] - mean[j]);

    /* fill the symmetric elements of the matrix. */
    mat[CMC(j, update, ncols)] = mat[CMC(update, j, ncols)] =
      (double)(temp / (nrows - 1));

  }/*FOR*/

}/*C_UPDATE_COVMAT*/

/* compute the sum of the squared errors. */
double c_sse(double *data, double mean, int nrows) {

int i = 0;
long double nvar = 0;

  for (i = 0; i < nrows; i++)
    nvar += (data[i] - mean) * (data[i] - mean);

  return (double)nvar;

}/*C_SSE*/

double c_mean(double *data, int nrows) {

int i = 0;
double mean = 0;

  for (i = 0; i < nrows; i++)
   mean += data[i];
  mean /= nrows;

  return mean;

}/*C_MEAN*/

void c_meanvec(double **data, double *mean, int nrows, int ncols, int first) {

int i = 0, j = 0;

  for (i = first; i < ncols; i++) {

    for (j = 0 ; j < nrows; j++)
      mean[i] += data[i][j];

    mean[i] /= nrows;

  }/*FOR*/

}/*C_MEANVEC*/

void c_ssevec(double **data, double *sse, double *means, int nrows, int ncols,
    int first) {

int i = 0, j = 0;
long double sum = 0;

  for (i = first; i < ncols; i++) {

    for (sum = 0, j = 0 ; j < nrows; j++)
      sum += (data[i][j] - means[i]) * (data[i][j] - means[i]);

    sse[i] = sum;

  }/*FOR*/

}/*C_SSEVEC*/

void c_update_meanvec(double **data, double *mean, int update, int nrows) {

int j = 0;

  mean[update] = 0;
  for (j = 0; j < nrows; j++)
    mean[update] += data[update][j];
  mean[update] /= nrows;

}/*C_UPDATE_MEANVEC*/
