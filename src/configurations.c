#include "include/rcore.h"
#include "include/sets.h"

/* wrapper around the c_configurations() function for use in .Call(). */
SEXP configurations(SEXP parents, SEXP factor, SEXP all) {

  return c_configurations(parents, isTRUE(factor), isTRUE(all));

}/*CFG2*/

/* wrapper around the cfg() function for use in C code. */
SEXP c_configurations(SEXP parents, int factor, int all_levels) {

int i = 0, *res = NULL, nlevels = 0;
SEXP temp, result;

  /* compute the configurations. */
  PROTECT(temp = allocVector(INTSXP, length(VECTOR_ELT(parents, 0))));
  cfg(parents, INTEGER(temp), &nlevels);

  if (factor) {

    /* convert the configurations from an integer array to to a factor. */
    if (all_levels)
      result = int2fac(temp, &nlevels);
    else
      result = int2fac(temp, NULL);

  }/*THEN*/
  else {

    /* keep the configurations as they are, but add 1 to get indexing right
     * in R. */
    result = temp;
    res = INTEGER(result);

    for (i = 0; i < length(result); i++)
      res[i]++;

  }/*ELSE*/

  UNPROTECT(1);

  return result;

}/*C_CONFIGURATIONS*/

/* identify different configurations of factors and assign them unique integer
 * codes, computed using the general formula for the column-mayor indexing. */
void cfg(SEXP parents, int *configurations, int *nlevels) {

int i = 0, **columns = NULL, *levels = NULL;
int ncols = length(parents), nrows = length(VECTOR_ELT(parents, 0));
SEXP temp;

  /* dereference the columns of the data frame. */
  columns = (int **) Calloc1D(ncols, sizeof(int *));
  levels = Calloc1D(ncols, sizeof(int));
  for (i = 0; i < ncols; i++) {

    temp = VECTOR_ELT(parents, i);
    columns[i] = INTEGER(temp);
    levels[i] = NLEVELS(temp);

  }/*FOR*/

  c_fast_config(columns, nrows, ncols, levels, configurations, nlevels, 0);

  Free1D(columns);
  Free1D(levels);

}/*CFG*/

void c_fast_config(int **columns, int nrows, int ncols, int *levels, int *configurations,
    int *nlevels, int offset) {

int i = 0, j = 0, cfgmap = 0;
long long *cumlevels = NULL, nl = 0;

  /* create the cumulative products of the number of levels. */
  cumlevels = Calloc1D(ncols, sizeof(long long));

  /* set the first one to 1 ... */
  cumlevels[0] = 1;

  /* ... then compute the following ones. */
  for (j = 1; j < ncols; j++)
    cumlevels[j] = cumlevels[j - 1] * levels[j - 1];

  /* compute the number of possible configurations. */
  nl = cumlevels[ncols - 1] * levels[ncols - 1];

  if (nl >= INT_MAX)
    error("attempting to create a factor with more than INT_MAX levels.");

  /* if nlevels is not a NULL pointer, save the number of possible
   * configurations. */
  if (nlevels)
    *nlevels = nl;

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
  configurations[i] = cfgmap + offset;

  }/*FOR*/

  Free1D(cumlevels);

}/*C_FAST_CONFIG*/

