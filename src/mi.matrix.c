#include "common.h"

#define DISCRETE_MAXIMUM_LIKELIHOOD 1
#define GAUSSIAN_MAXIMUM_LIKELIHOOD 2
#define DISCRETE_NETWORK(x) (x < GAUSSIAN_MAXIMUM_LIKELIHOOD)

#define LIST_MUTUAL_INFORMATION_COEFS() \
  if (*debuglevel > 0) { \
    for (i = 0; i < ncols; i++) \
      for (j = i + 1; j < ncols; j++) \
        Rprintf("  > mutual information between %s and %s is %lf.\n", \
          NODE(i), NODE(j), mim[UPTRI3(i + 1, j + 1, ncols)]); \
  }

#define CONVERT_TO_ARC_SET(ind, drop, rows) \
  PROTECT(arcs = allocMatrix(STRSXP, (rows), 2)); \
  for (i = 0, k = 0; i < ncols; i++) { \
    for (j = i + 1; j < ncols; j++) { \
      if (ind[UPTRI3(i + 1, j + 1, ncols)] != drop) { \
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, i)); \
         SET_STRING_ELT(arcs, k + (rows), STRING_ELT(nodes, j)); \
         k++; \
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, j)); \
         SET_STRING_ELT(arcs, k + (rows), STRING_ELT(nodes, i)); \
         k++; \
      } \
    } \
  } \
  finalize_arcs(arcs); \
  UNPROTECT(1);

#define DEREFERENCE_DATA_FRAME() \
  if (DISCRETE_NETWORK(*est)) { \
    columns = alloc1dpointer(ncols); \
    for (i = 0; i < ncols; i++) \
      columns[i] = INTEGER(VECTOR_ELT(data, i)); \
    nlevels = alloc1dcont(ncols); \
    for (i = 0; i < ncols; i++) \
      nlevels[i] = NLEVELS2(data, i); \
  } \
  else { \
    columns = (void **) alloc1dpointer(ncols); \
    for (i = 0; i < ncols; i++) \
      columns[i] = REAL(VECTOR_ELT(data, i)); \
  }

/* compute all the pairwise mutual information coefficients between the variables. */
void mi_matrix(double *mim, void **columns, int dim, int *nlevels, int *num,
    void *cond, int *clevels, int *est) {

int i = 0, j = 0;

  switch(*est) {

    case DISCRETE_MAXIMUM_LIKELIHOOD:

      if (cond == NULL) {

        for (i = 0; i < dim; i++) {

          for (j = i + 1; j < dim; j++) {

            mim[UPTRI3(i + 1, j + 1, dim)] =
              c_mi(((int **)columns)[i], nlevels + i,
                   ((int **)columns)[j], nlevels + j, num);

          }/*FOR*/

        }/*FOR*/

      }/*THEN*/
      else {

        for (i = 0; i < dim; i++) {

          for (j = i + 1; j < dim; j++) {

            mim[UPTRI3(i + 1, j + 1, dim)] =
              c_cmi(((int **)columns)[i], nlevels + i,
                    ((int **)columns)[j], nlevels + j,
                    (int *)cond, clevels, num);

          }/*FOR*/

        }/*FOR*/

      }/*ELSE*/

      break;

    case GAUSSIAN_MAXIMUM_LIKELIHOOD:

      for (i = 0; i < dim; i++) {

        for (j = i + 1; j < dim; j++) {

          mim[UPTRI3(i + 1, j + 1, dim)] =
            c_mig(((double **)columns)[i], ((double **)columns)[j], num);

        }/*FOR*/

      }/*FOR*/

    break;

  }/*SWITCH*/

}/*MI_MATRIX*/

/* ARACNE structure learning algorithm. */
SEXP aracne(SEXP data, SEXP estimator, SEXP whitelist, SEXP blacklist, SEXP debug) {

int i = 0, j = 0, k = 0, coord = 0, ncols = LENGTH(data);
int num = LENGTH(VECTOR_ELT(data, i)), narcs = ncols * (ncols - 1) / 2;
int *nlevels = NULL, *est = INTEGER(estimator), *wl = NULL, *bl = NULL;
int *debuglevel = LOGICAL(debug);
void **columns = NULL;
short int *exclude = NULL;
double *mim = NULL;
SEXP arcs, nodes, wlist, blist;

  nodes = getAttrib(data, R_NamesSymbol);

  /* dereference the columns of the data frame. */
  DEREFERENCE_DATA_FRAME()

  /* allocate the mutual information matrix and the status vector. */
  mim = alloc1dreal(UPTRI3_MATRIX(ncols));
  exclude = allocstatus(UPTRI3_MATRIX(ncols));

  /* compute the pairwise mutual information coefficients. */
  if (*debuglevel > 0)
    Rprintf("* computing pairwise mutual information coefficients.\n");

  mi_matrix(mim, columns, ncols, nlevels, &num, NULL, NULL, est);

  LIST_MUTUAL_INFORMATION_COEFS()

  /* compare all the triplets. */
  for (i = 0; i < ncols; i++) {

    for (j = i + 1; j < ncols; j++) {

      for (k = 0; k < ncols; k++) {

        if ((k == i) || (k == j))
          continue;

        /* cache the UPTRI3 coordinates of the arc. */
        coord = UPTRI3(i + 1, j + 1, ncols);

        /* if MI(X, Y) < min(MI(X, Z), MI(Z, Y)) drop arc X - Y. */
        if ((mim[coord] < mim[UPTRI3(i + 1, k + 1, ncols)]) &&
            (mim[coord] < mim[UPTRI3(j + 1, k + 1, ncols)])) {

          if (*debuglevel > 0) {

            Rprintf("* dropping arc %s - %s because of %s, %lf < min(%lf, %lf)\n",
              NODE(i), NODE(j), NODE(k), mim[UPTRI3(i + 1, j + 1, ncols)],
              mim[UPTRI3(i + 1, k + 1, ncols)], mim[UPTRI3(j + 1, k + 1, ncols)]);

          }/*THEN*/

          /* update the status vector. */
          exclude[coord] = 1;
          /* decrement the number of arcs. */
          narcs--;

          break;

        }/*THEN*/

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  /* add back whitelisted arcs. */
  if ((!isNull(whitelist)) && (LENGTH(whitelist) > 0)) {

    PROTECT(wlist = arc_hash(whitelist, nodes, TRUE, TRUE));
    wl = INTEGER(wlist);

    for (i = 0; i < LENGTH(wlist); i++) {

      if (*debuglevel > 0) {

        Rprintf("* adding back whitelisted arcs.\n");

        if (exclude[wl[i]] == 1) {

          Rprintf("  > arc %s - %s has been added to the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + LENGTH(wlist))));

        }/*THEN*/
        else {

          Rprintf("  > arc %s - %s was already present in the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + LENGTH(wlist))));

        }/*ELSE*/

      }/*THEN*/

      /* update the counter if need be. */
      if (exclude[wl[i]] == 1)
        narcs++;
      /* include the arc in the graph. */
      exclude[wl[i]] = 0;

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/

  /* remove blacklisted arcs. */
  if ((!isNull(blacklist)) && (LENGTH(blacklist) > 0)) {

    PROTECT(blist = arc_hash(blacklist, nodes, TRUE, TRUE));
    bl = INTEGER(blist);

    for (i = 0; i < LENGTH(blist); i++) {

      if (*debuglevel > 0) {

        Rprintf("* removing blacklisted arcs.\n");

        if (exclude[bl[i]] == 0) {

          Rprintf("  > arc %s - %s has been dropped from the graph.\n",
            CHAR(STRING_ELT(blacklist, i)), CHAR(STRING_ELT(blacklist, i + LENGTH(blist))));

        }/*THEN*/
        else {

          Rprintf("  > arc %s - %s was not present in the graph.\n",
            CHAR(STRING_ELT(blacklist, i)), CHAR(STRING_ELT(blacklist, i + LENGTH(blist))));

        }/*ELSE*/

      }/*THEN*/

      /* update the counter if need be. */
      if (exclude[bl[i]] == 0)
        narcs--;
      /* remove the arc from the graph. */
      exclude[bl[i]] = 1;

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/

  CONVERT_TO_ARC_SET(exclude, 1, 2 * narcs);

  return arcs;

}/*ARACNE*/

static int chow_liu_blacklist(int *blacklist, int *length, int *hash) {

  for (int k = 0; k < *length; k++)
    if (*hash == blacklist[k])
      return TRUE;

  return FALSE;

}/*CHOW_LIU_BLACKLIST*/

/* Chow-Liu structure learning algorithm. */
SEXP chow_liu(SEXP data, SEXP nodes, SEXP estimator, SEXP whitelist,
    SEXP blacklist, SEXP conditional, SEXP debug) {

int i = 0, j = 0, k = 0, debug_coord[2], ncols = LENGTH(data);
int num = LENGTH(VECTOR_ELT(data, i)), narcs = 0, nwl = 0, nbl = 0;
int *nlevels = NULL, *clevels = NULL, *est = INTEGER(estimator);
int *wl = NULL, *bl = NULL, *poset = NULL, *debuglevel = LOGICAL(debug);
void **columns = NULL, *cond = NULL;
short int *include = NULL;
double *mim = NULL;
SEXP arcs, wlist, blist;

  /* dereference the columns of the data frame. */
  DEREFERENCE_DATA_FRAME()

  /* only TAN uses a conditional variable, so assume it's discrete and go ahead. */
  if (conditional != R_NilValue) {

    cond = (void *) INTEGER(conditional);
    clevels = alloc1dcont(1);
    *clevels = NLEVELS(conditional); 

  }/*THEN*/

  /* allocate the mutual information matrix and the status vector. */
  mim = alloc1dreal(UPTRI3_MATRIX(ncols));
  include = allocstatus(UPTRI3_MATRIX(ncols));

  /* compute the pairwise mutual information coefficients. */
  if (*debuglevel > 0)
    Rprintf("* computing pairwise mutual information coefficients.\n");

  mi_matrix(mim, columns, ncols, nlevels, &num, cond, clevels, est);

  LIST_MUTUAL_INFORMATION_COEFS()

  /* add whitelisted arcs first. */
  if ((!isNull(whitelist)) && (LENGTH(whitelist) > 0)) {

    PROTECT(wlist = arc_hash(whitelist, nodes, TRUE, TRUE));
    wl = INTEGER(wlist);
    nwl = LENGTH(wlist);

    for (i = 0; i < nwl; i++) {

      if (*debuglevel > 0) {

        Rprintf("* adding whitelisted arcs first.\n");

        if (include[wl[i]] == 0) {

          Rprintf("  > arc %s - %s has been added to the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + nwl)));

        }/*THEN*/
        else {

          Rprintf("  > arc %s - %s was already present in the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + nwl)));

        }/*ELSE*/

      }/*THEN*/

      /* update the counter if need be. */
      if (include[wl[i]] == 0)
        narcs++;
      /* include the arc in the graph. */
      include[wl[i]] = 1;

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/

  /* cache blacklisted arcs. */
  if ((!isNull(blacklist)) && (LENGTH(blacklist) > 0)) {

    PROTECT(blist = arc_hash(blacklist, nodes, TRUE, TRUE));
    bl = INTEGER(blist);
    nbl = LENGTH(blist);

  }/*THEN*/

  /* sort the mutual information coefficients and keep track of the elements' index.  */
  poset = alloc1dcont(UPTRI3_MATRIX(ncols));
  for (i = 0; i < UPTRI3_MATRIX(ncols); i++)
    poset[i] = i;
  R_qsort_I(mim, poset, 1, UPTRI3_MATRIX(ncols));

  for (i = UPTRI3_MATRIX(ncols) - 1; i > 0; i--) {

    /* get back the coordinates from the position in the half-matrix. */
    INV_UPTRI3(poset[i], ncols, debug_coord);

    /* already included all the arcs we had to, exiting. */
    if (narcs >= ncols - 1)
      break;
    /* arc already present in the graph, nothing to do. */
    if (include[poset[i]] == 1)
      continue;

    if (bl) {

      if (chow_liu_blacklist(bl, &nbl, poset + i)) {

        if (*debuglevel > 0) {

          Rprintf("* arc %s - %s is blacklisted, skipping.\n",
            NODE(debug_coord[0]), NODE(debug_coord[1]));

        }/*THEN*/

        continue;

      }/*THEN*/

    }/*THEN*/

    if (c_uptri3_path(include, debug_coord[0], debug_coord[1], ncols, nodes, FALSE)) {

      if (*debuglevel > 0) {

        Rprintf("* arc %s - %s introduces cycles, skipping.\n",
          NODE(debug_coord[0]), NODE(debug_coord[1]));

      }/*THEN*/

      continue;

    }/*THEN*/

    if (*debuglevel > 0) {

      Rprintf("* adding arc %s - %s with mutual information %lf.\n",
        NODE(debug_coord[0]), NODE(debug_coord[1]), mim[i]);

    }/*THEN*/

    /* include the arc in the graph. */
    include[poset[i]] = 1;
    /* update the counter. */
    narcs++;

  }/*FOR*/

  if ((!isNull(blacklist)) && (LENGTH(blacklist) > 0))
    UNPROTECT(1);

  /* sanity check for blacklist-related madnes. */
  if (narcs != ncols - 1)
    error("learned %d arcs instead of %d, this is not a tree spanning all the nodes.",
      narcs, ncols - 1);

  CONVERT_TO_ARC_SET(include, 0, 2 * (ncols - 1));

  return arcs;

}/*CHOW_LIU*/

/* set the directions of the arcs in a tree given the root node. */
SEXP tree_directions(SEXP arcs, SEXP nodes, SEXP root, SEXP debug) {

int i = 0, j = 0, d = 0, traversed = 1;
int narcs = LENGTH(arcs)/2, nnodes = LENGTH(nodes);
int *a = NULL, *depth = 0, *debuglevel = LOGICAL(debug);
SEXP try, try2, result;

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  a = INTEGER(try);

  /* match the root node. */
  PROTECT(try2 = match(nodes, root, 0));

  /* allocate and initialize the statust vector. */
  depth = alloc1dcont(nnodes);
  depth[INT(try2) - 1] = 1;

  if (*debuglevel > 0)
    Rprintf("> root node (depth 1) is %s.\n", NODE(INT(try2) - 1));

  for (d = 1; d <= nnodes; d++) {

    if (*debuglevel > 0)
      Rprintf("> considering nodes at depth %d.\n", d + 1);

    for (i = 0; i < narcs; i++) {

      for (j = 0; j < nnodes; j++) {

        /* disregard nodes at the wrong depth. */
        if (depth[j] != d)
          continue;

        if ((a[i + narcs] == (j + 1)) && (depth[a[i] - 1] == 0)) {

          if (*debuglevel > 0)
            Rprintf("  * found node %s.\n", NODE(a[i] - 1));

          /* save the depth at which the node was found. */
          depth[a[i] - 1] = d + 1;

          /* update the counter of the traversed nodes. */
          traversed++;

        }/*THEN*/

      }/*FOR*/

    }/*FOR*/

    /* check whether all nodes have been traversed. */
    if (traversed == nnodes)
      break;

  }/*FOR*/

  /* allocate and initialize the return value. */
  PROTECT(result = allocMatrix(STRSXP, narcs/2, 2));

  for (i = 0, j = 0; i < narcs; i++) {

    if (depth[a[i] - 1] < depth[a[i + narcs] - 1]) {

      SET_STRING_ELT(result, j, STRING_ELT(arcs, i));
      SET_STRING_ELT(result, j + narcs/2, STRING_ELT(arcs, i + narcs));
      j++;

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(3);

  return result;

}/*TREE_DIRECTIONS*/
