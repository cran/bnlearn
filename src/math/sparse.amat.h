#ifndef SPARSE_AMAT_HEADER
#define SPARSE_AMAT_HEADER

#include <Rinternals.h>
#include <stdbool.h>

/* compressed sparse row representation of a square adjacency matrix. */
typedef struct {

  int dim;        /* number of rows, equal to the number of columns. */
  int narcs;      /* number of nonzero entries (directed arcs). */
  int *rowptr;    /* per-row offsets into colidx (length dim + 1). */
  int *colidx;    /* concatenated column indexes (length narcs). */
  char **labels;  /* row/column labels (length dim), or NULL. */

} sparse_amat;

sparse_amat new_sparse_amat(int *from, int *to, int narcs, int nnodes);
sparse_amat sparse_amat_from_SEXP(SEXP bn, bool transpose);
void FreeSparseAMAT(sparse_amat sam);

#endif
