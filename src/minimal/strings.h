#ifndef R_STRINGS_HEADER
#define R_STRINGS_HEADER

SEXP string_delete(SEXP array, SEXP string, int *idx);
SEXP string_setdiff(SEXP large, SEXP small);
SEXP mkStringVec(int n, ...);

#endif
