#ifndef UPPER_TRIANGULAR_HEADER
#define UPPER_TRIANGULAR_HEADER

/* upper triangular matrix (not including the diagonal). */
typedef struct {

  int dim;             /* dimension of the symmetrix matrix. */
  const char **names;  /* row and column names (assumed to be identical). */
  double *mat;         /* pointer to the matrix. */

} uppertriangular;

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

/* coordinate system for an upper triangular matrix (not including
 *  * the diagonal elements). */
#define UPTRI3(r, c, n) (UPTRI(r, c, n) - ((r > c) ? c : r))
void INV_UPTRI3(int x, int n, int *res);

/* dimension of the upper triangular part of a n x n matrix. */
#define UPTRI_MATRIX(n) (n) * ((n) + 1) / 2

/* upper triangular matrices. */
uppertriangular new_uppertriangular(int dim);
void uppertriangular_copy_names(uppertriangular *sym, const char **names);
void FreeUPPERTRIANGULAR(uppertriangular sym);
int uppertriangular_size(uppertriangular sym);
#define UTREL(sym, i, j) sym.mat[UPTRI3((i) + 1, (j) + 1, sym.dim)]

#endif
