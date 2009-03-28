#include <R.h>
#include <Rinternals.h>

/* numerical constants */

#define MACHINE_TOL 2.220446e-16

/* utility macros */

#define isTRUE(logical) LOGICAL(logical)[0] == TRUE
#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]

/* coordinate systems conversion matrices */

/*
 *  Coordinate system for an upper triangular matrix:
 *
 *  [(row - 1) * ncols + ncols] - [row * (row - 1) / 2]
 *
 *  the first term is the standard row major order coordinates;
 *  the second one is an adjustment to account for the missing
 *  lower half of the matrix.
 *
 */

/* this macro swaps its arguments to avoid "memory not mapped" errors. */
#define UPTRI(x, y, n) \
  (((x) <= (y)) ? \
    ((x) - 1) * n + (y) - ((x) * ((x) - 1)) / 2 : \
    ((y) - 1) * n + (x) - ((y) * ((y) - 1)) / 2)

/* this macro trusts its arguments to be correct, beware. */
#define UPTRI2(x, y, n) \
  ((x) - 1) * n + (y) - ((x) * ((x) - 1)) / 2

/* column-major coordinates for an arbitrary matrix. */
#define CMC(i, j, nrows) ((i) + (j) * (nrows))

/* utility functions  */

SEXP getListElement(SEXP list, char *str);

void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)

/* from arcs2amat.c */

SEXP arcs2amat(SEXP arcs, SEXP nodes);
SEXP amat2arcs(SEXP amat, SEXP nodes);

/* from cache.structure.c */

SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug);
SEXP cache_node_structure(int cur, SEXP nodes, SEXP amat, int nrows,
    int *status, SEXP debug);

/* from linear.algebra.c */

SEXP r_svd(SEXP matrix);

/* from linear.correlation.c */

SEXP fast_pcor(SEXP data, SEXP length);

/* from path.c */
SEXP c_has_dag_path(int start, int stop, int *amat, int n, SEXP nodes,
    SEXP underlying, SEXP debug);

/* memory allocation functions */

int *alloc1dcont(int length);
int **alloc2dcont(int length, int width);
int ***alloc3dcont(int length, int width, int depth);
short int *allocstatus(int length);
double *alloc1dreal(int length);
