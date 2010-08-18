#include "common.h"

/* populate the first subset (in lexicographic order). */
void first_subset(int *work, int *n) {

  for (int i = 0; i < *n; i++)
    work[i] = i;

}/*FIRST_SUBSET*/

/* get the next subset (in lexicographic order). */
int *next_subset(int *work, int *n, int *max) {

int i = 0, j = 0;

  /* this is the last possible subset, nothing to do. */
  if (work[0] == *max - *n)
    return NULL;

  if (work[*n - 1] < *max - 1) {

    /* increment the first element of the subset. */
    work[*n - 1]++;

  }/*THEN*/
  else {

    for (i = *n - 1; i >= 0; i--) {

      /* this is the last element of the set, look into the next slot... */
      if (work[i - 1] < work[i] - 1) {

        /* ... increment it... */
        work[i - 1]++;

        /* ... and reset the previous ones. */
        for (j = i - 1; j < *n; j++)
          work[j + 1] = work[j] + 1;

        break;

      }/*THEN*/

    }/*FOR*/

  }/*ELSE*/

  return work;

}/*NEXT_SUBSET*/

