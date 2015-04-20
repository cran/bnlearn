
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
 *  * the diagonal elements). */
#define UPTRI3(r, c, n) (UPTRI(r, c, n) - ((r > c) ? c : r))
void INV_UPTRI3(int x, int n, int *res);

/* dimension of the upper triangular part of a n x n matrix. */
#define UPTRI_MATRIX(n) (n) * ((n) + 1) / 2

/* dimension of the upper triangular part of a n x n matrix (not
 *  * including the diagonal elements). */
#define UPTRI3_MATRIX(n) (n) * ((n) - 1) / 2

/* column-major coordinates for an arbitrary matrix. */
#define CMC(i, j, nrows) ((i) + (j) * (nrows))

