#include "../include/rcore.h"

/* quick sort of an integer vector, with or without indexes. */
void i_sort(int *array, int *indexes, int length) {

  if (length == 0)
    return;

  if (!indexes)
    R_qsort_int(array, 1, length);
  else
    R_qsort_int_I(array, indexes, 1, length);

}/*I_SORT*/

/* quick sort of a double vector, with or without indexes. */
void d_sort(double *array, int *indexes, int length) {

  if (length == 0)
    return;

  if (!indexes)
    R_qsort(array, 1, length);
  else
    R_qsort_I(array, indexes, 1, length);

}/*D_SORT*/
