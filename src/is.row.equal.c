#include <R.h>
#include <Rinternals.h>

#define ARC(i,col) CHAR(STRING_ELT(data, i + col * nrows))

SEXP is_row_equal(SEXP data, SEXP array) {

  int i = 0;
  int nrows = LENGTH(data) / 2;
  const char *from = CHAR(STRING_ELT(array, 0));
  const char *to = CHAR(STRING_ELT(array, 1));
  SEXP res;

  PROTECT(res = allocVector(LGLSXP, nrows));

  for (i = 0; i < nrows; i++) {

    /* check the first element; if it does no match skip the second one. */
    if (!strcmp(from, ARC(i, 0)) ) {

      /* if the first element matches, check the other one. */
      if (!strcmp(to, ARC(i, 1)) ) {

        LOGICAL(res)[i] = TRUE;

      }/*THEN*/
      else {

        LOGICAL(res)[i] = FALSE;

      }/*ELSE*/

    }/*THEN*/
    else {

      LOGICAL(res)[i] = FALSE;

    }/*ELSE*/

  }/*FOR*/

  UNPROTECT(1);

  return res;

}/*IS_ROW_EQUAL*/

