#include "include/rcore.h"
#include "include/matrix.h"
#include "include/bn.h"

/* faster rbind() implementation for arc sets. */
SEXP arcs_rbind (SEXP matrix1, SEXP matrix2, SEXP reverse2) {

int i = 0, j = 0, m1 = length(matrix1)/2, m2 = length(matrix2)/2;
SEXP res;

  /* allocate the return value. */
  PROTECT(res = allocMatrix(STRSXP, m1 + m2, 2));
  /* allocate and initialize the column names. */
  setDimNames(res, R_NilValue, mkStringVec(2, "from", "to"));

  /* copy the elements of the first matrix. */
  for (i  = 0; i < m1; i++)
    for (j = 0; j < 2; j++)
      SET_STRING_ELT(res, CMC(i, j, m1 + m2), STRING_ELT(matrix1, CMC(i, j, m1)));

  /* copy the elements of the second matrix, reversing the order of the
   * columns as needed. */
  if (isTRUE(reverse2)) {

    for (i = 0; i < m2; i++)
      for(j = 0; j < 2; j++)
        SET_STRING_ELT(res, CMC(i + m1, j, m1 + m2), STRING_ELT(matrix2, CMC(i, 1 - j, m2)));

  }/*THEN*/
  else {

    for (i = 0; i < m2; i++)
      for(j = 0; j < 2; j++)
        SET_STRING_ELT(res, CMC(i + m1, j, m1 + m2), STRING_ELT(matrix2, CMC(i, j, m2)));

  }/*ELSE*/

  UNPROTECT(1);

  return res;

}/*ARCS_RBIND*/

