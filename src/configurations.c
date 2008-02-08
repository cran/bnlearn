#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/* .Call("cfg", data = alarm) */

/* macro for the [i,j] element of the data frame. */
#define BNLEARN_DATAFRAME(i, j) INTEGER(VECTOR_ELT(parents, j))[i]
/* macro for the number of levels of the [,j] column. */
#define BNLEARN_NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(parents, j), R_LevelsSymbol))

SEXP cfg(SEXP parents) {

  int i = 0, j = 0;
  int ncols = LENGTH(parents);
  int nrows = LENGTH(VECTOR_ELT(parents, 0));
  int cfgmap = 0;
  int cumlevels[nrows];
  
  SEXP ret;

  /* create the cumulative products of the number of levels. */
  memset(cumlevels, '\0', sizeof(int) * ncols);

  /* set the first one to 1 ... */
  cumlevels[0] = 1;

  /* ... then compute the following ones. */
  for (j = 1; j < ncols; j++) 
    cumlevels[j] = cumlevels[j - 1] * BNLEARN_NLEVELS(j - 1);

  /* allocate an array of size nrow for the configuration. */ 
  PROTECT(ret = allocVector(INTSXP, nrows));

  for (i = 0; i < nrows; i++) {

    /* reset the configuration mapping of the new row. */
    cfgmap = 0;

    for (j = 0; j < ncols; j++) {

      cfgmap += (BNLEARN_DATAFRAME(i, j) - 1) * cumlevels[j];

    }/*FOR*/

  /* save the configuration in the array. */
  INTEGER(ret)[i] = cfgmap;

  }/*FOR*/

  UNPROTECT(1);

return ret;

}/*CFG*/

