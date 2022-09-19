#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/uppertriangular.h"
#include "../../include/graph.h"
#include "../../tests/tests.h"
#include "../../include/bn.h"
#include "../../core/moments.h"
#include "../../core/correlation.h"
#include "../../core/data.table.h"
#include "../../core/contingency.tables.h"
#include "../../minimal/strings.h"
#include "../../minimal/common.h"

/* enum for the mutual information estimators, to be matched from the label
 * string passed down from R. */
typedef enum {
  ENOMI        =  0, /* error code, no such estimator. */

  MLE          =  1, /* maximum likelihood, discrete data. */

  MLE_G        = 10, /* maximum likelihood, Gaussian data. */
} mi_estimator_e;

#define ENTRY(key, value) if (strcmp(label, key) == 0) return value;

mi_estimator_e mi_to_enum(const char *label) {

  ENTRY("mi", MLE);
  ENTRY("mi-g", MLE_G);

  return ENOMI;

}/*MI_TO_ENUM*/

/* compute all the pairwise mutual information coefficients between discrete
 * variables. */
static void mi_matrix_discrete(uppertriangular mim, ddata data, int *cond,
    int clevels, mi_estimator_e est, bool debugging) {

int i = 0, j = 0;

  switch (est) {

    case MLE:

      if (!cond) {

        /* along the lines of Hartemink's discretisation code. */
        int max_nlvl = 0;
        counts2d counts = { 0 };

        /* figure out what is the largest contingency table that we need, so
         * we can allocate only one and shrink it as needed. */
        for (i = 0; i < data.m.ncols; i++)
          max_nlvl = (max_nlvl < data.nlvl[i]) ? data.nlvl[i] : max_nlvl;

        counts = new_2d_table(max_nlvl, max_nlvl, TRUE);

        for (i = 0; i < mim.dim; i++) {

          for (j = i + 1; j < mim.dim; j++) {

            /* resize the contigency table and fill it with the new counts. */
            resize_2d_table(data.nlvl[i], data.nlvl[j], &counts);
            refill_2d_table(data.col[i], data.col[j], &counts, data.m.nobs);

            /* compute and save the mutual information. */
            if (counts.nobs == 0)
              UTREL(mim, i, j) = 0;
            else
              UTREL(mim, i, j) = mi_kernel(counts) / counts.nobs;

            if (debugging)
              Rprintf("  > mutual information between %s and %s is %lf.\n",
                data.m.names[i], data.m.names[j], UTREL(mim, i, j));

          }/*FOR*/

        }/*FOR*/

        resize_2d_table(max_nlvl, max_nlvl, &counts);
        Free2DTAB(counts);

      }/*THEN*/
      else {

        /* along the lines of Hartemink's discretisation code. */
        int max_nlvl = 0;
        counts3d counts = { 0 };

        /* figure out what is the largest contingency table that we need, so
         * we can allocate only one and shrink it as needed. */
        for (i = 0; i < data.m.ncols; i++)
          max_nlvl = (max_nlvl < data.nlvl[i]) ? data.nlvl[i] : max_nlvl;

        counts = new_3d_table(max_nlvl, max_nlvl, clevels);

        for (i = 0; i < mim.dim; i++) {

          for (j = i + 1; j < mim.dim; j++) {

            /* resize the contigency table and fill it with the new counts. */
            resize_3d_table(data.nlvl[i], data.nlvl[j], clevels, &counts);
            refill_3d_table(data.col[i], data.col[j], cond, &counts, data.m.nobs);

            /* compute and save the mutual information. */

            if (counts.nobs == 0)
              UTREL(mim, i, j) = 0;
            else
              UTREL(mim, i, j) = cmi_kernel(counts) / counts.nobs;

            if (debugging)
              Rprintf("  > mutual information between %s and %s is %lf.\n",
                data.m.names[i], data.m.names[j], UTREL(mim, i, j));

          }/*FOR*/

        }/*FOR*/

        resize_3d_table(max_nlvl, max_nlvl, clevels, &counts);
        Free3DTAB(counts);

      }/*ELSE*/
      break;

    default:
      error("uknown mutual information estimator.");

  }/*SWITCH*/

}/*MI_MATRIX_DISCRETE*/

/* compute all the pairwise mutual information coefficients between Gaussian
 * variables. */
static void mi_matrix_gaussian(uppertriangular mim, gdata data, double *sse,
    mi_estimator_e est, bool debugging) {

int i = 0, j = 0;
double cor = 0;

  switch (est) {

    case MLE_G:

      for (i = 0; i < mim.dim; i++) {

        for (j = i + 1; j < mim.dim; j++) {

          if (data.m.flag[i].complete && data.m.flag[j].complete) {

            cor = c_fast_cor(data.col[i], data.col[j], data.m.nobs,
                    data.mean[i], data.mean[j], sse[i], sse[j]);

          }/*THEN*/
          else {

            cor = c_cor_with_missing(data.col[i], data.col[j], data.m.nobs,
                    NULL, NULL, NULL, NULL, NULL);

          }/*ELSE*/

          UTREL(mim, i, j) = cor_mi_trans(cor);

          if (debugging)
            Rprintf("  > mutual information between %s and %s is %lf.\n",
              data.m.names[i], data.m.names[j], UTREL(mim, i, j));

        }/*FOR*/

      }/*FOR*/
      break;

    default:
      error("uknown mutual information estimator.");

  }/*SWITCH*/

}/*MI_MATRIX_GAUSSIAN*/

uppertriangular estimate_mi_matrix(SEXP data, SEXP complete, SEXP conditional,
    mi_estimator_e est, bool debugging) {

int ncol = length(data);

uppertriangular mimatrix = { 0 };

  /* allocate the mutual information matrix. */
  mimatrix = new_uppertriangular(ncol);

  if (debugging)
    Rprintf("* computing pairwise mutual information coefficients.\n");

  if (est == MLE) {

    int clevels = 0, *cond = NULL;

    /* extract the data. */
    ddata dt = ddata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete, R_NilValue);
    /* the variable names are the dimension names of the mutual information
     * matrix. */
    uppertriangular_copy_names(&mimatrix, dt.m.names);

    /* only TAN uses a conditional variable, so assume it's discrete and go ahead. */
    if (conditional != R_NilValue) {

      cond = INTEGER(conditional);
      clevels = NLEVELS(conditional);

    }/*THEN*/

    /* compute the pairwise mutual information coefficients. */
    mi_matrix_discrete(mimatrix, dt, cond, clevels, est, debugging);

    FreeDDT(dt);

  }/*THEN*/
  else if (est == MLE_G) {

    /* extract the data, cache means and variances . */
    double *sse = NULL;
    gdata dt = gdata_from_SEXP(data, 0);
    meta_copy_names(&(dt.m), 0, data);
    meta_init_flags(&(dt.m), 0, complete, R_NilValue);
    gdata_cache_means(&dt, 0);
    sse = Calloc1D(ncol, sizeof(double));
    c_ssevec(dt.col, sse, dt.mean, dt.m.nobs, dt.m.ncols, 0);
    /* the variable names are the dimension names of the mutual information
     * matrix. */
    uppertriangular_copy_names(&mimatrix, dt.m.names);

    /* compute the pairwise mutual information coefficients. */
    mi_matrix_gaussian(mimatrix, dt, sse, est, debugging);

    Free1D(sse);
    FreeGDT(dt);

  }/*THEN*/
  else {

    error("uknown mutual information estimator.");

  }/*THEN*/

  return mimatrix;

}/*ESTIMATE_MI_MATRIX*/

/* ARACNE structure learning algorithm. */
SEXP aracne(SEXP data, SEXP estimator, SEXP whitelist, SEXP blacklist,
    SEXP complete, SEXP debug) {

int i = 0, j = 0, k = 0, ncol = length(data), narcs = ncol * (ncol - 1) / 2;
int *wl = NULL, *bl = NULL;
short int *exclude = NULL;
bool debugging = isTRUE(debug);
mi_estimator_e est = mi_to_enum(CHAR(STRING_ELT(estimator, 0)));
SEXP arcs, nodes, wlist, blist;
uppertriangular mimatrix = { 0 };

  PROTECT(nodes = getAttrib(data, R_NamesSymbol));

  /* estimate the mutual information matrix. */
  mimatrix = estimate_mi_matrix(data, complete, R_NilValue, est, debugging);

  exclude = Calloc1D(uppertriangular_size(mimatrix), sizeof(short int));

  /* compare all the triplets. */
  for (i = 0; i < ncol; i++) {

    for (j = i + 1; j < ncol; j++) {

      for (k = 0; k < ncol; k++) {

        if ((k == i) || (k == j))
          continue;

        /* if MI(X, Y) < min(MI(X, Z), MI(Z, Y)) drop arc X - Y. */
        if ((UTREL(mimatrix, i, j) < UTREL(mimatrix, i, k)) &&
            (UTREL(mimatrix, i, j) < UTREL(mimatrix, j, k))) {

          if (debugging) {

            Rprintf("* dropping arc %s - %s because of %s, %lf < min(%lf, %lf)\n",
              mimatrix.names[i], mimatrix.names[j], mimatrix.names[k],
              UTREL(mimatrix, i, j), UTREL(mimatrix, i, k), UTREL(mimatrix, j, k));

          }/*THEN*/

          /* update the status vector. */
          exclude[UPTRI3(i + 1, j + 1, ncol)] = 1;
          /* decrement the number of arcs. */
          narcs--;

          break;

        }/*THEN*/

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  /* add back whitelisted arcs. */
  if ((!isNull(whitelist)) && (length(whitelist) > 0)) {

    PROTECT(wlist = arc_hash(whitelist, nodes, TRUE, TRUE));
    wl = INTEGER(wlist);

    for (i = 0; i < length(wlist); i++) {

      if (debugging) {

        Rprintf("* adding back whitelisted arcs.\n");

        if (exclude[wl[i]] == 1) {

          Rprintf("  > arc %s - %s has been added to the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + length(wlist))));

        }/*THEN*/
        else {

          Rprintf("  > arc %s - %s was already present in the graph.\n",
            CHAR(STRING_ELT(whitelist, i)), CHAR(STRING_ELT(whitelist, i + length(wlist))));

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
  if ((!isNull(blacklist)) && (length(blacklist) > 0)) {

    PROTECT(blist = arc_hash(blacklist, nodes, TRUE, TRUE));
    bl = INTEGER(blist);

    for (i = 0; i < length(blist); i++) {

      if (debugging) {

        Rprintf("* removing blacklisted arcs.\n");

        if (exclude[bl[i]] == 0) {

          Rprintf("  > arc %s - %s has been dropped from the graph.\n",
            CHAR(STRING_ELT(blacklist, i)), CHAR(STRING_ELT(blacklist, i + length(blist))));

        }/*THEN*/
        else {

          Rprintf("  > arc %s - %s was not present in the graph.\n",
            CHAR(STRING_ELT(blacklist, i)), CHAR(STRING_ELT(blacklist, i + length(blist))));

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

  PROTECT(arcs = allocMatrix(STRSXP, (2 * narcs), 2));
  for (i = 0, k = 0; i < ncol; i++) {
    for (j = i + 1; j < ncol; j++) {
      if (exclude[UPTRI3(i + 1, j + 1, ncol)] != 1) {
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, i));
         SET_STRING_ELT(arcs, k + (2 * narcs), STRING_ELT(nodes, j));
         k++;
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, j));
         SET_STRING_ELT(arcs, k + (2 * narcs), STRING_ELT(nodes, i));
         k++;
      }/*THEN*/
    }/*FOR*/
  }/*FOR*/
  setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));
  UNPROTECT(1);

  FreeUPPERTRIANGULAR(mimatrix);
  Free1D(exclude);

  UNPROTECT(1);

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
    SEXP blacklist, SEXP complete, SEXP conditional, SEXP debug) {

int i = 0, j = 0, k = 0, debug_coord[2], ncol = length(data);
int narcs = 0, nwl = 0, nbl = 0, *depth = NULL;
int *wl = NULL, *bl = NULL, *poset = NULL;
short int *include = NULL;
bool debugging = isTRUE(debug);
mi_estimator_e est = mi_to_enum(CHAR(STRING_ELT(estimator, 0)));
SEXP arcs, wlist, blist;
uppertriangular mimatrix = { 0 };

  /* estimate the mutual information matrix. */
  mimatrix = estimate_mi_matrix(data, complete, conditional, est, debugging);

  include = Calloc1D(uppertriangular_size(mimatrix), sizeof(short int));

  /* add whitelisted arcs first. */
  if ((!isNull(whitelist)) && (length(whitelist) > 0)) {

    PROTECT(wlist = arc_hash(whitelist, nodes, TRUE, TRUE));
    wl = INTEGER(wlist);
    nwl = length(wlist);

    for (i = 0; i < nwl; i++) {

      if (debugging) {

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
  if ((!isNull(blacklist)) && (length(blacklist) > 0)) {

    PROTECT(blist = arc_hash(blacklist, nodes, TRUE, TRUE));
    bl = INTEGER(blist);
    nbl = length(blist);

  }/*THEN*/

  /* sort the mutual information coefficients and keep track of the elements' index.  */
  poset = Calloc1D(uppertriangular_size(mimatrix), sizeof(int));
  for (i = 0; i < uppertriangular_size(mimatrix); i++)
    poset[i] = i;
  R_qsort_I(mimatrix.mat, poset, 1, uppertriangular_size(mimatrix));

  depth = Calloc1D(ncol, sizeof(int));

  for (i = uppertriangular_size(mimatrix) - 1; i >= 0; i--) {

    /* get back the coordinates from the position in the half-matrix. */
    INV_UPTRI3(poset[i], ncol, debug_coord);

    /* already included all the arcs we had to, exiting. */
    if (narcs >= ncol - 1)
      break;
    /* arc already present in the graph, nothing to do. */
    if (include[poset[i]] == 1)
      continue;

    if (bl) {

      if (chow_liu_blacklist(bl, &nbl, poset + i)) {

        if (debugging) {

          Rprintf("* arc %s - %s is blacklisted, skipping.\n",
            NODE(debug_coord[0]), NODE(debug_coord[1]));

        }/*THEN*/

        continue;

      }/*THEN*/

    }/*THEN*/

    if (c_uptri3_path(include, depth, debug_coord[0], debug_coord[1], ncol,
          nodes, FALSE)) {

      if (debugging) {

        Rprintf("* arc %s - %s introduces cycles, skipping.\n",
          NODE(debug_coord[0]), NODE(debug_coord[1]));

      }/*THEN*/

      continue;

    }/*THEN*/

    if (debugging) {

      Rprintf("* adding arc %s - %s with mutual information %lf.\n",
        NODE(debug_coord[0]), NODE(debug_coord[1]), mimatrix.mat[i]);

    }/*THEN*/

    /* include the arc in the graph. */
    include[poset[i]] = 1;
    /* update the counter. */
    narcs++;

  }/*FOR*/

  if ((!isNull(blacklist)) && (length(blacklist) > 0))
    UNPROTECT(1);

  /* sanity check for blacklist-related madnes. */
  if (narcs != ncol - 1)
    error("learned %d arcs instead of %d, this is not a tree spanning all the nodes.",
      narcs, ncol - 1);

  PROTECT(arcs = allocMatrix(STRSXP, (2 * (ncol - 1)), 2));
  for (i = 0, k = 0; i < ncol; i++) {
    for (j = i + 1; j < ncol; j++) {
      if (include[UPTRI3(i + 1, j + 1, ncol)] != 0) {
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, i));
         SET_STRING_ELT(arcs, k + (2 * (ncol - 1)), STRING_ELT(nodes, j));
         k++;
         SET_STRING_ELT(arcs, k, STRING_ELT(nodes, j));
         SET_STRING_ELT(arcs, k + (2 * (ncol - 1)), STRING_ELT(nodes, i));
         k++;
      }/*THEN*/
    }/*FOR*/
  }/*FOR*/
  setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));
  UNPROTECT(1);

  Free1D(depth);
  FreeUPPERTRIANGULAR(mimatrix);
  Free1D(include);
  Free1D(poset);

  return arcs;

}/*CHOW_LIU*/

/* set the directions of the arcs in a tree given the root node. */
SEXP tree_directions(SEXP arcs, SEXP nodes, SEXP root, SEXP debug) {

int i = 0, j = 0, d = 0, traversed = 1;
int narcs = length(arcs)/2, nnodes = length(nodes);
int *a = NULL, *depth = 0;
bool debugging = isTRUE(debug);
SEXP try, try2, result;

  /* match the node labels in the arc set. */
  PROTECT(try = match(nodes, arcs, 0));
  a = INTEGER(try);

  /* match the root node. */
  PROTECT(try2 = match(nodes, root, 0));

  /* allocate and initialize the statust vector. */
  depth = Calloc1D(nnodes, sizeof(int));
  depth[INT(try2) - 1] = 1;

  if (debugging)
    Rprintf("> root node (depth 1) is %s.\n", NODE(INT(try2) - 1));

  for (d = 1; d <= nnodes; d++) {

    if (debugging)
      Rprintf("> considering nodes at depth %d.\n", d + 1);

    for (i = 0; i < narcs; i++) {

      for (j = 0; j < nnodes; j++) {

        /* disregard nodes at the wrong depth. */
        if (depth[j] != d)
          continue;

        if ((a[i + narcs] == (j + 1)) && (depth[a[i] - 1] == 0)) {

          if (debugging)
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

  Free1D(depth);

  return result;

}/*TREE_DIRECTIONS*/
