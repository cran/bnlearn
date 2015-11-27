#include "include/rcore.h"
#include "include/globals.h"

/* a rudimental C implementation of which.max() for a double array. */
int d_which_max(double *array, int length) {

int i = 0, imax = -1;
double max = R_NegInf;

  for (i = 0; i < length; i++) {

    /* NA and NaN cannot be compared with valid real numbers. */
    if (ISNAN(array[i]))
      continue;

    if (array[i] > max) {

      imax = i;
      max = array[i];

    }/*THEN*/

  }/*FOR*/

  if (imax < 0) {

    /* if all elements are -Inf, return the first as the maximum. */
    if (array[0] == R_NegInf)
      return 1;

    /* if all elements are NA/NaN return NA. */
    return NA_INTEGER;

  }/*THEN*/

  return imax + 1;

}/*D_WHICH_MAX*/

/* a rudimental C implementation of which.max() for a long double array. */
int ld_which_max(long double *array, int length) {

int i = 0, imax = -1;
long double max = R_NegInf;

  for (i = 0; i < length; i++) {

    /* NA and NaN cannot be compared with valid real numbers. */
    if (ISNAN(array[i]))
      continue;

    if (array[i] > max) {

      imax = i;
      max = array[i];

    }/*THEN*/

  }/*FOR*/

  if (imax < 0) {

    /* if all elements are -Inf, return the first as the maximum. */
    if (array[0] == R_NegInf)
      return 1;

    /* if all elements are NA/NaN return NA. */
    return NA_INTEGER;

  }/*THEN*/

  return imax + 1;

}/*LD_WHICH_MAX*/

/* return all maxima in the array, modulo numeric tolerance. */
int all_max(double *array, int length, int *maxima, int *indexes,
    double *buf) {

int i = 0, nmax = 0;
double tol = MACHINE_TOL;

  /* make a safety copy of the array. */
  memcpy(buf, array, length * sizeof(double));

  /* sort the elements of the array. */
  rsort_with_index(buf, indexes, length);

  /* count the number of maxima (considering numeric tolerance). */
  for (i = length - 1; i >= 0; i--)
    if (buf[i] < buf[length - 1] - tol)
      break;

  /* set the counter for the number of maxima. */
  nmax = length - i - 1;

  /* save the indexes of the maxima. */
  memcpy(maxima, indexes + length - nmax, nmax * sizeof(int));

  return nmax;

}/*ALL_MAX*/

