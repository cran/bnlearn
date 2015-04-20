#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Arith.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>

/* utility macros */
#define isTRUE(logical) LOGICAL(logical)[0] == TRUE
#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]
#define NODE(i) CHAR(STRING_ELT(nodes, i))

/* macro for the number of levels of the [,j] column. */
#define NLEVELS(x) length(getAttrib(x, R_LevelsSymbol))
#define NLEVELS2(data, j) \
  length(getAttrib(VECTOR_ELT(data, j), R_LevelsSymbol))

/* from common.c */
SEXP getListElement(SEXP list, char *str);
void *DATAPTR(SEXP x);
SEXP unique(SEXP array);
SEXP dupe(SEXP array);
SEXP mkRealVec(int n, ...);
SEXP int2fac(SEXP vector, int *nlevels);

/* from strings.c */
SEXP string_delete(SEXP array, SEXP string, int *idx);
SEXP mkStringVec(int n, ...);

/* from which.max.c */
int all_max(double *array, int length, int *maxima, int *indexes, double *buf);
int d_which_max(double *array, int length);
int ld_which_max(long double *array, int length);

