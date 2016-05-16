#include "include/rcore.h"
#include "include/matrix.h"

/* populate the first subset (in lexicographic order). */
void first_subset(int *work, int n, int offset) {

  for (int i = 0; i < n; i++)
    work[i] = i + offset;

}/*FIRST_SUBSET*/

/* get the next subset (in lexicographic order). */
int next_subset(int *work, int n, int max, int offset) {

int i = 0, j = 0;

  /* this is the last possible subset, nothing to do. */
  if (work[0] - offset == max - n)
    return FALSE;

  if (work[n - 1] - offset < max - 1) {

    /* increment the first element of the subset. */
    work[n - 1]++;

  }/*THEN*/
  else {

    for (i = n - 1; i >= 0; i--) {

      /* this is the last element of the set, look into the next slot... */
      if (work[i - 1] < work[i] - 1) {

        /* ... increment it... */
        work[i - 1]++;

        /* ... and reset the previous ones. */
        for (j = i - 1; j < n - 1; j++)
          work[j + 1] = work[j] + 1;

        break;

      }/*THEN*/

    }/*FOR*/

  }/*ELSE*/

  return TRUE;

}/*NEXT_SUBSET*/

/* enumerate all subsets of a certain size (R interface). */
SEXP r_subsets(SEXP elems, SEXP size) {

int i = 0, k = 0, n = length(elems), r = INT(size), *id = NULL;
double nsub = choose(n, r);
SEXP result;

 if (nsub * r > INT_MAX)
   error("too many subsets of size %d.", r);

 /* allocate the scratch space and the return value. */
 id = Calloc1D(r, sizeof(int));
 PROTECT(result = allocMatrix(STRSXP, nsub, r));

  /* iterate over subsets. */
  first_subset(id, r, 0);

  for (k = 0;  k < nsub; k++) {

    for (i = 0; i < r; i++)
      SET_STRING_ELT(result, CMC(k, i, nsub), STRING_ELT(elems, id[i]));

    next_subset(id, r, n, 0);

  }/*FOR*/

  Free1D(id);
  UNPROTECT(1);

  return result;

}/*R_SUBSETS*/
