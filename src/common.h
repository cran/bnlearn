#include <R.h>
#include <Rinternals.h>

/* numerical constants */

#define MACHINE_TOL sqrt(DOUBLE_EPS)

/* utility macros */

#define isTRUE(logical) LOGICAL(logical)[0] == TRUE
#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]
#define NODE(i) CHAR(STRING_ELT(nodes, i))

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
    ((x) - 1) * n + (y) - 1 - ((x) * ((x) - 1)) / 2 : \
    ((y) - 1) * n + (x) - 1 - ((y) * ((y) - 1)) / 2)

/* this macro trusts its arguments to be correct, beware. */
#define UPTRI2(x, y, n) \
  ((x) - 1) * n + (y) - 1 - ((x) * ((x) - 1)) / 2

/* coordinate system for an upper triangular matrix (not including
 * the diagonal elements). */
#define UPTRI3(r, c, n) (UPTRI(r, c, n) - ((r > c) ? c : r))

/* dimension of the upper triangular part of a n x n matrix. */
#define UPTRI_MATRIX(n) (n) * ((n) + 1) / 2

/* dimension of the upper triangular part of a n x n matrix (not
 * including the diagonal elements). */
#define UPTRI3_MATRIX(n) (n) * ((n) - 1) / 2

/* column-major coordinates for an arbitrary matrix. */
#define CMC(i, j, nrows) ((i) + (j) * (nrows))

/* utility functions  */

SEXP getListElement(SEXP list, char *str);
SEXP unique(SEXP array);
SEXP dupe(SEXP array);
int which_max(double *array, int length);
SEXP finalize_arcs(SEXP arcs); 

void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)

SEXP int2fac(SEXP vector);

/* from arcs2amat.c */

SEXP arcs2amat(SEXP arcs, SEXP nodes);
SEXP amat2arcs(SEXP amat, SEXP nodes);

/* from cache.structure.c */

SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug);
SEXP c_cache_partial_structure(int target, SEXP nodes, SEXP amat, int *status, SEXP debug);

/* from linear.algebra.c */

SEXP r_svd(SEXP matrix);
SEXP r_det(SEXP matrix, int scale);
double c_det(double *matrix, int *rows);

/* from linear.correlation.c */

SEXP fast_pcor(SEXP data, SEXP length);

/* from path.c */

int c_has_path(int start, int stop, int *amat, int n, SEXP nodes,
    int ugraph, int notdirect, int debuglevel);
int c_directed_path(int start, int stop, int *amat, int n, SEXP nodes,
    int debuglevel);

/* from hash.c */

SEXP arc_hash(SEXP arcs, SEXP nodes);
SEXP c_amat_hash(int *amat, int *nnodes);

/* from configurations.c */

SEXP cfg(SEXP parents);

/* shared between hill climbing and tabu search. */

void bestop_update(SEXP bestop, char *op, const char *from, const char *to);

/* from filter.arcs.c */

SEXP c_unique_arcs(SEXP arcs, SEXP nodes, int warnlevel);

/* memory allocation functions */

int *alloc1dcont(int length);
int **alloc2dcont(int length, int width);
int ***alloc3dcont(int length, int width, int depth);
short int *allocstatus(int length);
double *alloc1dreal(int length);
double **alloc2dreal(int length, int width);
double ***alloc3dreal(int length, int width, int depth);
char **alloc1dstring (int length);
