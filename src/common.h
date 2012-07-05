#include <R.h>
#include <Rinternals.h>

/* numerical constants */

#define MACHINE_TOL sqrt(DOUBLE_EPS)

/* utility macros */

#define isTRUE(logical) LOGICAL(logical)[0] == TRUE
#define INT(x) INTEGER(x)[0]
#define NUM(x) REAL(x)[0]
#define NODE(i) CHAR(STRING_ELT(nodes, i))

/* macro for the number of levels of the [,j] column. */
#define NLEVELS(x) LENGTH(getAttrib(x, R_LevelsSymbol))
#define NLEVELS2(data, j) \
  LENGTH(getAttrib(VECTOR_ELT(data, j), R_LevelsSymbol))

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
void INV_UPTRI3(int x, int n, int *res);

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
void all_max(double *array, int length, int *maxima, int *nmax, int *indexes);
SEXP finalize_arcs(SEXP arcs);
SEXP minimal_data_frame(SEXP obj);
SEXP dataframe_column(SEXP dataframe, SEXP name, SEXP drop);
SEXP int2fac(SEXP vector, int *nlevels);

/* from sampling.c */

void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)
void SampleReplace(int k, int n, int *y, int *x);
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans);

/* from arcs2amat.c */

SEXP arcs2amat(SEXP arcs, SEXP nodes);
SEXP amat2arcs(SEXP amat, SEXP nodes);

/* from cache.structure.c */

SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug);
SEXP c_cache_partial_structure(int target, SEXP nodes, SEXP amat, int *status, SEXP debug);

/* from linear.algebra.c */

SEXP r_svd(SEXP matrix, SEXP strict);
SEXP r_det(SEXP matrix, int scale);
double c_det(double *matrix, int *rows);
void c_svd(double *A, double *U, double *D, double *V, int *nrows, int *ncols,
    int *mindim, int strict, int *errcode);
void c_ginv(double *covariance, int *ncols, double *mpinv);
void c_finv(double *cov, int *ncols, double *mpinv);
double c_quadratic(double *x, int *ncols, double *sigma, double *y,
    double *workspace);
void c_rotate(double *S1, double *S2, double *x, double *a, double *mu,
    int *ncols, double *workspace);

/* from linear.correlation.c */

void c_covmat(double **data, double *mean, int *ncols, int *nrows, double *mat);
void c_update_covmat(double **data, double *mean, int update, int *ncols,
    int *nrows, double *mat);
double c_fast_cor(double *xx, double *yy, int *num);
double c_fast_pcor(double *covariance, int *ncols, double *u, double *d,
    double *vt, int *errcode);
SEXP fast_pcor(SEXP data, SEXP length, SEXP shrinkage, SEXP strict);

/* from shrinkage.c */

SEXP cov_lambda(SEXP data, SEXP length);

/* from path.c */

int c_has_path(int start, int stop, int *amat, int n, SEXP nodes,
    int ugraph, int notdirect, int debuglevel);
int c_directed_path(int start, int stop, int *amat, int n, SEXP nodes,
    int debuglevel);
int c_uptri3_path(short int *uptri, int from, int to, int nnodes,
    SEXP nodes, int debuglevel);

/* from hash.c */

SEXP arc_hash(SEXP arcs, SEXP nodes, int uptri, int sort);
SEXP c_amat_hash(int *amat, int *nnodes);

/* from configurations.c */

void cfg(SEXP parents, int *configurations, int *nlevels);

/* shared between hill climbing and tabu search. */

void bestop_update(SEXP bestop, char *op, const char *from, const char *to);

/* from filter.arcs.c */

SEXP c_unique_arcs(SEXP arcs, SEXP nodes, int warnlevel);

/* from fitted.c */

SEXP root_nodes(SEXP bn, SEXP leaves);

/* from simulation.c */

SEXP schedule(SEXP bn, SEXP root_nodes, SEXP reverse, SEXP debug);

/* from mutual.information.c */

#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))

double c_mi(int *xx, int *llx, int *yy, int *lly, int *num);
double c_cmi(int *xx, int *llx, int *yy, int *lly, int *zz, int *llz, int *num);
double c_mig(double *xx, double *yy, int *num);

/* memory allocation functions */

int *alloc1dcont(int length);
int **alloc2dcont(int length, int width);
int ***alloc3dcont(int length, int width, int depth);
short int *allocstatus(int length);
double *alloc1dreal(int length);
double **alloc2dreal(int length, int width);
double ***alloc3dreal(int length, int width, int depth);
void **alloc1dpointer (int length);
char **alloc1dstring (int length);
