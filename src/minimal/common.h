#ifndef R_COMMON_HEADER
#define R_COMMON_HEADER

/* macro for the number of levels of a factor. */
#define NLEVELS(x) length(getAttrib(x, R_LevelsSymbol))

/* from common.c */
bool c_is(SEXP obj, const char *str);
SEXP getListElement(SEXP list, char *str);
void *DATAPTR(SEXP x);
SEXP mkReal(double x);
SEXP mkRealVec(int n, ...);
SEXP int2fac(SEXP vector, int *nlevels);
void setDimNames(SEXP obj, SEXP rownames, SEXP colnames);
SEXP subset_by_name(SEXP vec, int n, ...);
bool all_equal(SEXP vec, SEXP val);

#endif
