#include "common.h"

/* wrapper around the cfg() function for use in .Call(). */
SEXP cfg2(SEXP parents, SEXP factor, SEXP all) {

int i = 0, *res = NULL, nlevels = 0, *all_levels = LOGICAL(all);
SEXP temp, result;

  /* compute the configurations. */
  PROTECT(temp = allocVector(INTSXP, LENGTH(VECTOR_ELT(parents, 0))));
  cfg(parents, INTEGER(temp), &nlevels);

  if (isTRUE(factor)) {

    /* convert the configurations from an integer array to to a factor. */
    if (*all_levels)
      result = int2fac(temp, &nlevels);
    else
      result = int2fac(temp, NULL);

  }/*THEN*/
  else {

    /* keep the configurations as they are, but add 1 to get indexing right
     * in R. */
    result = temp;
    res = INTEGER(result);

    for (i = 0; i < LENGTH(result); i++)
      res[i]++;

  }/*ELSE*/

  UNPROTECT(1);

  return result;

}/*CFG2*/

/* identify different configurations of factors and assign them unique integer
 * codes, computed using the general formula for the column-mayor indexing. */
void cfg(SEXP parents, int *configurations, int *nlevels) {

int i = 0, j = 0, cfgmap = 0;
int ncols = LENGTH(parents), nrows = LENGTH(VECTOR_ELT(parents, 0));
long long *cumlevels = NULL, nl = 0;
int **columns = NULL;

  /* create the cumulative products of the number of levels. */
  cumlevels = (long long *) R_alloc(ncols, sizeof(long long));
  memset(cumlevels, '\0', sizeof(long long) * ncols); 

  /* dereference the columns of the data frame. */
  columns = (int **) alloc1dpointer(ncols);
  for (i = 0; i < ncols; i++)
    columns[i] = INTEGER(VECTOR_ELT(parents, i));

  /* set the first one to 1 ... */
  cumlevels[0] = 1;

  /* ... then compute the following ones. */
  for (j = 1; j < ncols; j++)
    cumlevels[j] = cumlevels[j - 1] * NLEVELS2(parents, j - 1);

  /* compute the number of possible configurations. */
  nl = cumlevels[ncols - 1] * NLEVELS2(parents, ncols - 1);

  if (nl >= INT_MAX) {

    error("attempting to create a factor with more than INT_MAX levels.");

  }/*THEN*/
  else {

    /* if nlevels is not a NULL pointer, save the number of possible
      * configurations. */
    if (nlevels)
      *nlevels = nl;

  }/*ELSE*/

  for (i = 0; i < nrows; i++) {

    /* reset the configuration mapping of the new row. */
    cfgmap = 0;

    for (j = 0; j < ncols; j++) {

      if (columns[j][i] == NA_INTEGER) {

        cfgmap = NA_INTEGER;
        break;

      }/*THEN*/
      else {

        cfgmap += (columns[j][i] - 1) * cumlevels[j];

     }/*ELSE*/

    }/*FOR*/

  /* save the configuration in the array. */
  configurations[i] = cfgmap;

  }/*FOR*/

}/*CFG*/

