#include "include/rcore.h"
#include "include/globals.h"
#include "include/matrix.h"
#include "include/dataframe.h"
#include "include/graph.h"

SEXP smart_network_averaging(SEXP arcs, SEXP nodes, SEXP weights) {

int k = 0, from = 0, to = 0, nrows = length(arcs) / 2, dims = length(nodes);
int *a = NULL, *coords = NULL, *poset = NULL, *path = NULL, *scratch = NULL;
double *w = NULL;
SEXP weights2, amat, try, acyclic;

  /* allocate and initialize the adjacency matrix. */
  PROTECT(amat = allocMatrix(INTSXP, dims, dims));
  a = INTEGER(amat);
  memset(a, '\0', sizeof(int) * dims * dims);

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  coords = INTEGER(try);

  /* duplicate the weights to preserve the original ones. */
  PROTECT(weights2 = duplicate(weights));
  w = REAL(weights2);

  /* sort the strength coefficients. */
  poset = Calloc1D(nrows, sizeof(int));
  for (k = 0; k < nrows; k++)
    poset[k] = k;
  R_qsort_I(w, poset, 1, nrows);

  /* allocate buffers for c_has_path(). */
  path = Calloc1D(dims, sizeof(int));
  scratch = Calloc1D(dims, sizeof(int));

  /* iterate over the arcs in reverse order wrt their strength coefficients. */
  for (k = 0; k < nrows; k++) {

    from = coords[poset[k]] - 1;
    to = coords[poset[k] + nrows] - 1;

    /* add an arc only if it does not introduce cycles. */
    if (!c_has_path(to, from, a, dims, nodes, FALSE, TRUE, path, scratch, FALSE))
      a[CMC(from, to, dims)] = 1;
    else
      warning("arc %s -> %s would introduce cycles in the graph, ignoring.",
        NODE(from), NODE(to));

  }/*FOR*/

  /* convert the adjacency matrix back to an arc set and return it. */
  acyclic = amat2arcs(amat, nodes);

  Free1D(path);
  Free1D(scratch);
  Free1D(poset);

  UNPROTECT(3);

  return acyclic;

}/*SMART_NETWORK_AVERAGING*/

/* mean strength for p-values and score deltas. */
static void mean_strength_overall(SEXP *mean_df, SEXP strength, SEXP nodes,
  int nrows, int nstr, SEXP ref_hash, double *w) {

int i = 0, j = 0, *t = NULL;
double *mstr = NULL, *cur_strength = NULL;
long double cumw = 0;
SEXP mean_str, cur, cur_hash, try;

  /* allocate the strength accumulator vector. */
  PROTECT(mean_str = allocVector(REALSXP, nrows));
  SET_VECTOR_ELT(*mean_df, 2, mean_str);
  mstr = REAL(mean_str);
  memset(mstr, '\0', nrows * sizeof(double));

  for (i = 0; i < nstr; i++) {

    /* move to the next object. */
    cur = VECTOR_ELT(strength, i);
    /* get the strength values from the bn.strength object. */
    cur_strength = REAL(VECTOR_ELT(cur, 2));
    /* get the arc IDs to use to correctly match strengths. */
    PROTECT(cur_hash = arc_hash(cur, nodes, FALSE, FALSE));

    /* match the current arc IDs to the reference arc IDs. */
    PROTECT(try = match(ref_hash, cur_hash, 0));
    t = INTEGER(try);

    for (j = 0; j < nrows; j++)
      mstr[t[j] - 1] += w[i] * cur_strength[j];

    /* update the total weight mass. */
    cumw += w[i];

    UNPROTECT(2);

  }/*FOR*/

  /* rescale by the total weight mass. */
  for (j = 0; j < nrows; j++)
    mstr[j] /= cumw;

  UNPROTECT(1);

}/*MEAN_STRENGTH_OVERALL*/

/* mean strength for bootstrap probabilities. */
static void mean_strength_direction(SEXP *mean_df, SEXP strength, SEXP nodes,
  int nrows, int nstr, SEXP ref_hash, double *w) {

int i = 0, j = 0, *t = NULL, nnodes = length(nodes);
double *mstr = NULL, *mdir = NULL, *cur_strength = NULL, *cur_dir = NULL;
double fwd = 0, bkwd = 0;
long double cumw = 0;
SEXP mean_str, mean_dir, cur, cur_hash, try;


  /* allocate vectors for strength and direction. */
  PROTECT(mean_str = allocVector(REALSXP, nrows));
  SET_VECTOR_ELT(*mean_df, 2, mean_str);
  mstr = REAL(mean_str);
  memset(mstr, '\0', nrows * sizeof(double));
  PROTECT(mean_dir = allocVector(REALSXP, nrows));
  SET_VECTOR_ELT(*mean_df, 3, mean_dir);
  mdir = REAL(mean_dir);
  memset(mdir, '\0', nrows * sizeof(double));

  for (i = 0; i < nstr; i++) {

    /* move to the next object. */
    cur = VECTOR_ELT(strength, i);
    /* get the strength and direction values from the bn.strength object. */
    cur_strength = REAL(VECTOR_ELT(cur, 2));
    cur_dir = REAL(VECTOR_ELT(cur, 3));
    /* get the arc IDs to use to correctly match strengths. */
    PROTECT(cur_hash = arc_hash(cur, nodes, FALSE, FALSE));

    /* match the current arc IDs to the reference arc IDs. */
    PROTECT(try = match(ref_hash, cur_hash, 0));
    t = INTEGER(try);

    for (j = 0; j < nrows; j++)
      mstr[t[j] - 1] += w[i] * (cur_strength[j] * cur_dir[j]);

    /* update the total weight mass. */
    cumw += w[i];

    UNPROTECT(2);

  }/*FOR*/

  /* rescale by the total weight mass. */
  for (j = 0; j < nrows; j++)
    mstr[j] /= cumw;

  /* split arc strength from direction strength. */
  for (i = 0; i < nnodes; i++) {

    for (j = i + 1; j < nnodes; j++) {

      fwd = mstr[CMC(j, i, nnodes) - i - 1];
      bkwd = mstr[CMC(i, j, nnodes) - j];

      mstr[CMC(j, i, nnodes) - i - 1] = mstr[CMC(i, j, nnodes) - j] = fwd + bkwd;

      if (bkwd + fwd > 0) {

      mdir[CMC(j, i, nnodes) - i - 1] = fwd / (fwd + bkwd);
      mdir[CMC(i, j, nnodes) - j] = bkwd / (fwd + bkwd);

      }/*THEN*/
      else {

        mdir[CMC(j, i, nnodes) - i - 1] = mdir[CMC(i, j, nnodes) - j] = 0;

      }/*ELSE*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(2);

}/*MEAN_STRENGTH_DIRECTION*/

/* average multiple bn.strength objects, with weights. */
SEXP mean_strength(SEXP strength, SEXP nodes, SEXP weights) {

int nstr = length(weights), ncols = 0, nrows = 0;
double *w = REAL(weights);
const char *m = NULL;
SEXP ref, ref_hash, mean_df, method;

  /* initialize the result using the first bn.strength object as a reference. */
  ref = VECTOR_ELT(strength, 0);
  ncols = length(ref);
  nrows = length(VECTOR_ELT(ref, 0));

  PROTECT(mean_df = allocVector(VECSXP, ncols));
  setAttrib(mean_df, R_NamesSymbol, getAttrib(ref, R_NamesSymbol));
  SET_VECTOR_ELT(mean_df, 0, VECTOR_ELT(ref, 0));
  SET_VECTOR_ELT(mean_df, 1, VECTOR_ELT(ref, 1));
  /* make it a data frame */
  minimal_data_frame(mean_df);

  /* compute the arc IDs to match arcs of later bn.strength objects. */
  PROTECT(ref_hash = arc_hash(ref, nodes, FALSE, FALSE));

  /* switch backend according to how the strengths were computed. */
  method = getAttrib(ref, BN_MethodSymbol);
  m = CHAR(STRING_ELT(method, 0));

  if ((strcmp(m, "score") == 0) || (strcmp(m, "test") == 0))
    mean_strength_overall(&mean_df, strength, nodes, nrows, nstr, ref_hash, w);
  else if (strcmp(m, "bootstrap") == 0)
    mean_strength_direction(&mean_df, strength, nodes, nrows, nstr, ref_hash, w);

  UNPROTECT(2);

  return mean_df;

}/*MEAN_STRENGTH*/
