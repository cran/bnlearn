#include "../include/rcore.h"

void i_sort(int *array, int *indexes, int length) {

  if (length == 0)
    return;

  if (!indexes)
    R_qsort_int(array, 1, length);
  else
    R_qsort_int_I(array, indexes, 1, length);

}/*I_SORT*/

void d_sort(double *array, int *indexes, int length) {

  if (length == 0)
    return;

  if (!indexes)
    R_qsort(array, 1, length);
  else
    R_qsort_I(array, indexes, 1, length);

}/*D_SORT*/
