#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Arith.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>
#include <stdbool.h>

/* for backwards compatibility with older R versions. */
#ifndef MAYBE_REFERENCED
#define MAYBE_REFERENCED(x) (NAMED(x) > 0)
#endif

/* utility macros. */
#define isTRUE(logical) ((LOGICAL(logical)[0]) == TRUE)
#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]
#define NODE(i) CHAR(STRING_ELT(nodes, i))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

