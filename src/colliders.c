#include "include/rcore.h"
#include "include/graph.h"
#include "include/matrix.h"

/* return the colliders present in the graph. */
int c_colliders(int *a, int nnodes, int **colliders, bool want_shielded,
    bool want_unshielded, char **node_labels, bool debugging) {

int i = 0, j = 0, k = 0, cursor = 0, capacity = nnodes * 3;
bool shielded = FALSE;

  /* for each potential head of a directed arc... */
  for (j = 0; j < nnodes; j++) {

    if (debugging)
      Rprintf("* looking at arcs pointing at node %s.\n", node_labels[j]);

    /* ... look for a tail... */
    for (i = 0; i < nnodes; i++) {

      /* no directed arc here. */
      if (!((a[CMC(i, j, nnodes)] > 0) && (a[CMC(j, i, nnodes)] == 0)))
        continue;

      if (debugging)
        Rprintf("  > found arc %s -> %s.\n", node_labels[i], node_labels[j]);

      /* ... look for a second directed arc with the same head and a
       * different tail... */
      for (k = i + 1; k < nnodes; k++) {

        /* no directed arc here. */
        if (!((a[CMC(k, j, nnodes)] > 0) && (a[CMC(j, k, nnodes)] == 0)))
          continue;

        if (debugging)
          Rprintf("    > found a second arc %s -> %s.\n", node_labels[k],
            node_labels[j]);

        /* ... check whether the two tails are adjacent or not... */
        shielded = ((a[CMC(i, k, nnodes)] > 0) || (a[CMC(k, i, nnodes)] > 0));

        /* .. and save it if it is a collider we want to return. */
        if ((shielded && want_shielded) || (!shielded && want_unshielded)) {

          /* extend as needed. */
          if (cursor + 3 > capacity)
            Realloc1D(*colliders, capacity + nnodes * 3, sizeof(int));

          (*colliders)[cursor++] = i;
          (*colliders)[cursor++] = j;
          (*colliders)[cursor++] = k;

          if (debugging)
            Rprintf("    @ found %s collider %s -> %s <- %s.\n",
              (shielded) ? "shielded" : "unshielded",
              node_labels[i], node_labels[j], node_labels[k]);

        }/*THEN*/

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  return cursor / 3;

}/*C_COLLIDERS*/

/* R interface to c_colliders(). */
SEXP colliders(SEXP arcs, SEXP nodes, SEXP return_arcs, SEXP shielded,
    SEXP unshielded, SEXP debug) {

int i = 0, nnodes = length(nodes), ncoll = 0;
int *a = NULL, *coll = NULL;
char **node_labels = NULL;
SEXP amat, result;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* extract the node labels from the adjacency matrix. */
  node_labels = Calloc1D(nnodes, sizeof(char *));
  for (i = 0; i < nnodes; i++)
    node_labels[i] = (char *)CHAR(STRING_ELT(nodes, i));

  /* allocate the colliders as a row-major matrix with a default length,
   * to make it easy to make it bigger if needed. */
  coll = Calloc1D(nnodes * 3, sizeof(int));

  /* identify the colliders. */
  ncoll = c_colliders(a, nnodes, &coll, isTRUE(shielded), isTRUE(unshielded),
            node_labels, isTRUE(debug));

  /* contruct the return value. */
  PROTECT(result = allocMatrix(STRSXP, ncoll, 3));
  setDimNames(result, R_NilValue, mkStringVec(3, "X", "Z", "Y"));

  for (i = 0; i < ncoll; i++) {

    SET_STRING_ELT(result, CMC(i, 0, ncoll), STRING_ELT(nodes, coll[i * 3 + 0]));
    SET_STRING_ELT(result, CMC(i, 1, ncoll), STRING_ELT(nodes, coll[i * 3 + 1]));
    SET_STRING_ELT(result, CMC(i, 2, ncoll), STRING_ELT(nodes, coll[i * 3 + 2]));

  }/*FOR*/

  Free1D(coll);
  Free1D(node_labels);

  UNPROTECT(2);

  return result;

}/*COLLIDERS*/

