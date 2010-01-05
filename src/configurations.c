#include "common.h"

/* macro for the [i,j] element of the data frame. */
#define DATAFRAME(i, j) INTEGER(VECTOR_ELT(parents, j))[i]
/* macro for the number of levels of the [,j] column. */
#define NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(parents, j), R_LevelsSymbol))

SEXP cfg(SEXP parents);

SEXP cfg2(SEXP parents) {

SEXP temp, result;

  PROTECT(temp = cfg(parents));
  result = int2fac(temp);

  UNPROTECT(1);

  return result;

}/*CFG2*/

SEXP cfg(SEXP parents) {

int i = 0, j = 0, cfgmap = 0;
int ncols = LENGTH(parents), nrows = LENGTH(VECTOR_ELT(parents, 0));
int *cumlevels = NULL, *res = NULL;
SEXP result;

  /* create the cumulative products of the number of levels. */
  cumlevels = alloc1dcont(ncols);

  /* set the first one to 1 ... */
  cumlevels[0] = 1;

  /* ... then compute the following ones. */
  for (j = 1; j < ncols; j++)
    cumlevels[j] = cumlevels[j - 1] * NLEVELS(j - 1);

  /* allocate an array of size nrow for the configuration. */
  PROTECT(result = allocVector(INTSXP, nrows));
  res = INTEGER(result);

  for (i = 0; i < nrows; i++) {

    /* reset the configuration mapping of the new row. */
    cfgmap = 0;

    for (j = 0; j < ncols; j++) {

      cfgmap += (DATAFRAME(i, j) - 1) * cumlevels[j];

    }/*FOR*/

  /* save the configuration in the array. */
  res[i] = cfgmap;

  }/*FOR*/

  UNPROTECT(1);

return result;

}/*CFG*/

