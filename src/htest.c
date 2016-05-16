#include "include/rcore.h"
#include "include/tests.h"

SEXP c_create_htest(double stat, SEXP test, double pvalue, double df, SEXP B) {

test_e test_type = test_label(CHAR(STRING_ELT(test, 0)));
SEXP result, s, n, params;

  /* allocate the return value. */
  PROTECT(result = allocVector(VECSXP, 7));
  /* set the class. */
  setAttrib(result, R_ClassSymbol, mkString("htest"));
  /* set the names of the elements. */
  setAttrib(result, R_NamesSymbol, mkStringVec(7, "statistic", "p.value",
    "method", "null.value", "alternative", "data.name", "parameter"));

  /* set the test statistic. */
  PROTECT(s = ScalarReal(stat));
  setAttrib(s, R_NamesSymbol, test);
  SET_VECTOR_ELT(result, 0, s);

  /* set the p-value. */
  SET_VECTOR_ELT(result, 1, ScalarReal(pvalue));

  /* set the label of the test. */
  SET_VECTOR_ELT(result, 2, mkString(""));

  /* set the value of the statistic under the null. */
  PROTECT(n = ScalarReal(0));
  setAttrib(n, R_NamesSymbol, mkString("value"));
  SET_VECTOR_ELT(result, 3, n);

  /* set the alternative hypothesis. */
  if (IS_TWO_SIDED(test_type))
    SET_VECTOR_ELT(result, 4, mkString("two.sided"));
  else
    SET_VECTOR_ELT(result, 4, mkString("greater"));

  /* set the data description string. */
  SET_VECTOR_ELT(result, 5, mkString(""));

  if (ISNAN(df)) {

    if (B != R_NilValue) {

      PROTECT(params = ScalarReal(INT(B)));
      setAttrib(params, R_NamesSymbol, mkString("Monte Carlo samples"));
      SET_VECTOR_ELT(result, 6, params);
      UNPROTECT(1);

    }/*THEN*/

  }/*THEN*/
  else {

    if (B != R_NilValue) {

      PROTECT(params = allocVector(REALSXP, 2));
      REAL(params)[0] = df;
      REAL(params)[1] = INT(B);
      setAttrib(params, R_NamesSymbol, mkStringVec(2, "df", "Monte Carlo samples"));
      SET_VECTOR_ELT(result, 6, params);
      UNPROTECT(1);

    }/*THEN*/
    else {

      PROTECT(params = ScalarReal(df));
      setAttrib(params, R_NamesSymbol, mkString("df"));
      SET_VECTOR_ELT(result, 6, params);
      UNPROTECT(1);

    }/*ELSE*/

  }/*ELSE*/

  UNPROTECT(3);

  return result;

}/*C_CREATE_HTEST*/

