#include "include/rcore.h"
#include "include/graph.h"
#include "include/matrix.h"

#define ABSENT     0
#define PRESENT    1
#define FIXED      2
#define IMMUTABLE  3

#define UNCHANGED  0
#define CHANGED    1

static void scan_graph(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel);
static void mark_vstructures(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel);
static void unmark_shielded(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel);
static int prevent_cycles(int *a, SEXP nodes, int nnodes, int debuglevel);
static int prevent_additional_vstructures(int *a, SEXP nodes, int nnodes,
    int debuglevel);
static int prevent_conflicts(int *a, SEXP nodes, int nnodes, short int *scratch,
    int debuglevel);
static int prevent_chains(int *a, SEXP nodes, int nnodes, int debuglevel);
static void renormalize_amat(int *a, int nnodes);
static SEXP amat2vstructs(int *a, SEXP nodes, int nnodes, short int *collider);

/* return the v-structures present in the graph. */
SEXP vstructures(SEXP arcs, SEXP nodes, SEXP return_arcs, SEXP including_moral,
    SEXP debug) {

int i = 0, nnodes = length(nodes);
int *a = NULL, *ret_arcs = LOGICAL(return_arcs), debuglevel = isTRUE(debug);
short int *collider = NULL;
SEXP amat, result;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = Calloc1D(nnodes, sizeof(short int));

  /* STEP 1: scan the graph. */
  scan_graph(a, nodes, nnodes, collider, debuglevel);

  /* STEP 2: all the arcs not part of a v-structure are now undirected. */
  mark_vstructures(a, nodes, nnodes, collider, debuglevel);

  /* decide whether or not to keep moral v-structures (i.e. shielded colliders). */
  if (!isTRUE(including_moral))
    unmark_shielded(a, nodes, nnodes, collider, debuglevel);

  /* remove all non-fixed arcs (including compelled ones). */
  for (i = 0; i < nnodes * nnodes; i++)
    a[i] = (a[i] == FIXED) ? 1 : 0;

  if (*ret_arcs > 0)
    PROTECT(result = amat2arcs(amat, nodes));
  else
    PROTECT(result = amat2vstructs(a, nodes, nnodes, collider));

  UNPROTECT(2);

  Free1D(collider);

  return result;

}/*VSTRUCTURES*/

static SEXP amat2vstructs(int *a, SEXP nodes, int nnodes, short int *collider) {

int i = 0, j = 0, k = 0, nvstructs = 0, counter = 0, row = 0;
SEXP result;

  for (i = 0; i < nnodes; i++) {

    /* this is not a collider, nothing to do. */
    if (collider[i] == 0)
      continue;

    /* reset the parents' counter. */
    counter = 0;

    /* count the parents. */
    for (j = 0; j < nnodes; j++)
      counter += a[CMC(j, i, nnodes)];

    /* only arcs that are part of a v-structure are still present in the
     * adjacency matrix, so computing the number of v-structures is easy. */
    nvstructs += counter * (counter - 1) / 2;

  }/*FOR*/

  PROTECT(result = allocMatrix(STRSXP, nvstructs, 3));

  /* allocate and the colnames. */
  setDimNames(result, R_NilValue, mkStringVec(3, "X", "Z", "Y"));

  for (i = 0; i < nnodes; i++) {

    /* this is not a collider, nothing to do. */
    if (collider[i] == 0)
      continue;

    for (j = 0; j < nnodes; j++) {

      /* this is no parent. */
      if (a[CMC(j, i, nnodes)] == 0)
        continue;

      /* this is a parent, look for the other one. */
      for (k = j + 1; k < nnodes; k++) {

        /* again this is no parent. */
        if (a[CMC(k, i, nnodes)] == 0)
          continue;

        SET_STRING_ELT(result, CMC(row, 0, nvstructs), STRING_ELT(nodes, j));
        SET_STRING_ELT(result, CMC(row, 1, nvstructs), STRING_ELT(nodes, i));
        SET_STRING_ELT(result, CMC(row, 2, nvstructs), STRING_ELT(nodes, k));
        row++;

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(1);

  return result;

}/*AMAT2VSTRUCTS*/

static int fix_wlbl(int *a, int nnodes, SEXP nodes, SEXP whitelist,
  SEXP blacklist, int debuglevel) {

int i = 0, j = 0, *w = NULL, *b = NULL;
SEXP wl, bl;

  if ((whitelist == R_NilValue) && (blacklist == R_NilValue)) {

    /* skipping STEP 2, just propagate the directions from the arcs that are
     * directed (this is for learning algorithms, to preserve whitelisted and
     * blacklisted arcs). */

    for (i = 0; i < nnodes; i++) {

      for (j = 0; j < nnodes; j++) {

        if ((a[CMC(i, j, nnodes)] == PRESENT) && (a[CMC(j, i, nnodes)] == ABSENT))
          a[CMC(i, j, nnodes)] = FIXED;
        if ((a[CMC(i, j, nnodes)] == ABSENT) && (a[CMC(j, i, nnodes)] == PRESENT))
          a[CMC(j, i, nnodes)] = FIXED;

      }/*FOR*/

    }/*FOR*/

    return TRUE;

  }/*THEN*/
  else {

    /* do not skip STEP 2, but fix the directions in the whitelist and remove
       the arcs in the blacklist. */

    if (whitelist != R_NilValue) {

      PROTECT(wl = arcs2amat(whitelist, nodes));
      w = INTEGER(wl);

      for (i = 0; i < nnodes; i++) {
        for (j = 0; j < nnodes; j++) {

          /* if a directed arc is in the whitelist, fix its direction. */
          if ((a[CMC(i, j, nnodes)] == PRESENT) && (a[CMC(j, i, nnodes)] == ABSENT) &&
              (w[CMC(i, j, nnodes)] == PRESENT)) {

            if (debuglevel > 0)
              Rprintf("  > marking arc %s -> %s as immutable.\n", NODE(i), NODE(j));

            a[CMC(i, j, nnodes)] = IMMUTABLE;

          }/*THEN*/

        }/*FOR*/
      }/*FOR*/

      UNPROTECT(1);

    }/*THEN*/

    if (blacklist != R_NilValue) {

      PROTECT(bl = arcs2amat(blacklist, nodes));
      b = INTEGER(bl);

      for (i = 0; i < nnodes; i++) {
        for (j = 0; j < nnodes; j++) {

          /* if a directed arc is in the blacklist, and the opposite is in
           * the graph, fix the direction of the arc in the graph  */
          if ((a[CMC(i, j, nnodes)] == PRESENT) && (b[CMC(j, i, nnodes)] == PRESENT) &&
              (b[CMC(i, j, nnodes)] == ABSENT)) {

            if (debuglevel > 0)
              Rprintf("  > marking arc %s -> %s as immutable.\n", NODE(i), NODE(j));

            a[CMC(i, j, nnodes)] = IMMUTABLE;

          }/*THEN*/

        }/*FOR*/
      }/*FOR*/

      UNPROTECT(1);

    }/*THEN*/

  }/*ELSE*/

  return FALSE;

}/*FIX_WLBL*/

SEXP cpdag(SEXP arcs, SEXP nodes, SEXP moral, SEXP fix, SEXP whitelist,
    SEXP blacklist, SEXP illegal, SEXP debug) {

int i = 0, changed = 0, nnodes = length(nodes);
short int *collider = NULL, *scratch = NULL;
int *a = NULL, debuglevel = isTRUE(debug), all_fixed = FALSE;
SEXP amat;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = Calloc1D(nnodes, sizeof(short int));
  /* allocate a second scratch vector. */
  scratch = Calloc1D(nnodes, sizeof(short int));

  /* STEP 1a: scan the graph. */
  scan_graph(a, nodes, nnodes, collider, debuglevel);

  /* STEP 1b: fix arcs whose direction is determined by parametric assumptions
   * (the illegal argument works like a blacklist). */
  if (illegal != R_NilValue) {

   if (debuglevel > 0)
     Rprintf("* setting the directions determined by parametric assumptions (step 1b).\n");

    fix_wlbl(a, nnodes, nodes, R_NilValue, illegal, debuglevel);

  }/*THEN*/

  /* STEP 1c: fix arcs according to the whitelist and the blacklist. */
  if (isTRUE(fix)) {

   if (debuglevel > 0)
     Rprintf("* setting the directions determined by whitelist and blacklist (step 1c).\n");

    all_fixed = fix_wlbl(a, nnodes, nodes, whitelist, blacklist, debuglevel);

  }/*THEN*/

  if (!all_fixed) {

    /* STEP 2: all the arcs not part of a v-structure are now undirected. */
    mark_vstructures(a, nodes, nnodes, collider, debuglevel);

    /* decide whether or not to keep moral v-structures (i.e. shielded colliders). */
    if (isTRUE(moral))
      unmark_shielded(a, nodes, nnodes, collider, debuglevel);

  }/*THEN*/

  /* STEP 3: orient some more edges without introducing cycles or new
   * v-structures in the graph. */
  for (i = 0; i < nnodes * nnodes; i++) {

    if (debuglevel > 0)
      Rprintf("* setting the direction of more arcs (step 3, iteration %d).\n", i);

    /* STEP 3a: prevent the creation of cycles (Meek Rule 2). */
    changed |= prevent_cycles(a, nodes, nnodes, debuglevel);

    /* STEP 3b: prevent the creation of new v-structures (Meek Rule 1). */
    changed |= prevent_additional_vstructures(a, nodes, nnodes, debuglevel);

    /* STEP 3c: do not leave arcs undirected if they rule out some directions
     * of some other undirected arcs (Meek Rule 3). */
    changed |= prevent_conflicts(a, nodes, nnodes, scratch, debuglevel);

    /* STEP 3d: avoid being stuck due to long chains mixing directed and
     * undirected arcs (Meek Rule 4). */
    changed |= prevent_chains(a, nodes, nnodes, debuglevel);

    if (!changed) {

      if (debuglevel > 0)
        Rprintf("  > the graph is unchanged, stopping.\n");

      break;

    }/*THEN*/

    /* reset changed for the next iteration. */
    changed = 0;

  }/*FOR*/

  /* STEP 4: */

  /* renormalize the adjacency matrix (i.e. set all elements either to 0 or 1). */
  renormalize_amat(a, nnodes);

  UNPROTECT(1);

  Free1D(collider);
  Free1D(scratch);

  return amat;

}/*CPDAG*/

static void scan_graph(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel) {

 int i = 0, j = 0, counter = 0;

  /* count the parents of each node (the non-symmetric 1s in the
   * corresponding column of the adjacency matrix). */
  if (debuglevel > 0)
    Rprintf("* scanning the graph (step 1).\n");

  for (j = 0; j < nnodes; j++) {

    counter = 0;

    for (i = 0; i < nnodes; i++) {

      /* loops are never present in a (partially) directed acyclic graph. */
      if (i == j) continue;

      if ((a[CMC(i, j, nnodes)] == PRESENT) && (a[CMC(j, i, nnodes)] == ABSENT))
        counter++;

      /* this node has at least two parents, so it's the collider in a
       * v-structure. Mark it as such and skip to the next one. */
      if (counter >= 2) {

        collider[j] = 1;

        if (debuglevel > 0)
          Rprintf("  > node %s is the center of a v-structure.\n", NODE(j));

        break;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

}/*SCAN_GRAPH*/

static void mark_vstructures(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel) {

int i = 0, j = 0;

  /* set arcs belonging to a v-structure as FIXED. */
  if (debuglevel > 0)
    Rprintf("* marking v-structures (step 2).\n");

  for (j = 0; j < nnodes; j++) {

    /* no v-structure here, mark all (non IMMUTABLE) arcs as undirected and skip
     * to the next node. */
    if (!collider[j]) {

      for (i = 0; i < nnodes; i++) {

        if (a[CMC(i, j, nnodes)] == PRESENT)
          a[CMC(j, i, nnodes)] = PRESENT;

      }/*FOR*/

    }/*THEN*/

  }/*FOR*/

  for (j = 0; j < nnodes; j++) {

    /* no v-structure here, this node will be examined in step 3. */
    if (!collider[j]) continue;

    for (i = 0; i < nnodes; i++) {

      /* loops are not allowed in (partially) directed graphs. */
      if (i == j) continue;

      /* if the arc is marked as IMMUTABLE, it is always directed so nothing to do here. */
      if (a[CMC(i, j, nnodes)] == PRESENT) {

        if (a[CMC(j, i, nnodes)] == ABSENT) {

          /* this is a directed arc, it's part of a v-structure; mark it as FIXED. */
          a[CMC(i, j, nnodes)] = FIXED;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, it's part of a v-structure.\n", NODE(i), NODE(j));

        }/*THEN*/
        else if (a[CMC(j, i, nnodes)] == PRESENT) {

          if (collider[i]) {

            /* at this point the arcs is undirected, and each endpoint is the
             * center of a v-structure: this means that the arc was undirected
             * in the original graph, or it would have been fixed as a directed
             * arc in a previous step.*/
            a[CMC(i, j, nnodes)] = a[CMC(j, i, nnodes)] = FIXED;

            if (debuglevel > 0)
              Rprintf("  > fixing arc %s - %s due to conflicting v-structures.\n", NODE(j), NODE(i));

          }/*THEN*/

        }/*THEN*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

}/*MARK_VSTRUCTURES*/

static void unmark_shielded(int *a, SEXP nodes, int nnodes, short int *collider,
    int debuglevel) {

int i = 0, j = 0, k = 0, check = FALSE;

  if (debuglevel > 0)
    Rprintf("* removing moral v-structures aka shielded colliders.\n");

  for (i = 0; i < nnodes; i++) {

    /* that's not a collider, nothing to do.*/
    if (collider[i] == 0)
      continue;

    for (j = 0; j < nnodes; j++) {

      /* this arc is not part of a v-structure. */
      if ((a[CMC(j, i, nnodes)] != FIXED) && (a[CMC(j, i, nnodes)] != IMMUTABLE))
        continue;

      if (debuglevel > 0)
        Rprintf("  > considering arc %s -> %s.\n", NODE(j), NODE(i));

      /* this is set to TRUE is the arc is part of an unshielded collider. */
      check = FALSE;

      for (k = 0; k < nnodes; k++) {

        /* this arc is not part of a v-structure. */
        if (((a[CMC(k, i, nnodes)] != FIXED) && (a[CMC(k, i, nnodes)] != IMMUTABLE)) || (j == k))
          continue;

        if (debuglevel > 0)
          Rprintf("    > considering v-structure %s -> %s <- %s.\n", NODE(j), NODE(i), NODE(k));

        /* it's a real v-structure if the parents are not connected or if both
         * arcs are marked IMMUTABLE. */
        if (((a[CMC(j, k, nnodes)] == ABSENT) && (a[CMC(k, j, nnodes)] == ABSENT)) ||
            ((a[CMC(j, i, nnodes)] == IMMUTABLE) && (a[CMC(k, i, nnodes)] == IMMUTABLE))) {

          check = TRUE;

          break;

        }/*THEN*/

      }/*FOR*/

      if (check) {

        if (debuglevel > 0)
          Rprintf("  @ arc %s -> %s is part of a v-structure.\n", NODE(j), NODE(i));

      }/*THEN*/
      else {

        if (debuglevel > 0)
          Rprintf("  @ arc %s -> %s is not part of a v-structure.\n", NODE(j), NODE(i));

        /* this arc should not be marked as FIXED and should not be directed;
         * leave it alone if marked IMMUTABLE. */
        if (a[CMC(j, i, nnodes)] == FIXED)
          a[CMC(j, i, nnodes)] = a[CMC(i, j, nnodes)] = PRESENT;

      }/*ELSE*/

    }/*FOR*/

  }/*FOR*/

}/*UNMARK_SHIELDED*/

static int vstructures_and_cycles(int *a, int nnodes, int i, int j) {

int k = 0;

  for (k = 0; k < nnodes; k++) {

    /* no valid arcs. */
    if ((k == i) || (k == j))
      continue;

    /* looking for a v-structure centered in NODE(j), so for NODE(k) -> NODE(j)
     * and no arc between NODE(i) and NODE(j). */
    if ((a[CMC(k, j, nnodes)] >= PRESENT) && (a[CMC(j, k, nnodes)] == ABSENT) &&
        (a[CMC(k, i, nnodes)] == ABSENT) && (a[CMC(i, k, nnodes)] == ABSENT))
      return TRUE;

  }/*FOR*/

  return FALSE;

}/*VSTRUCTURES_VS_CYCLES*/

static int prevent_cycles(int *a, SEXP nodes, int nnodes, int debuglevel) {

int i = 0, j = 0, changed = UNCHANGED;
int *path = NULL, *scratch = NULL;

  /* allocate buffers for c_directed_path(). */
  path = Calloc1D(nnodes, sizeof(int));
  scratch = Calloc1D(nnodes, sizeof(int));

  for (j = 0; j < nnodes; j++) {

    for (i = 0; i < nnodes; i++) {

      if ((a[CMC(i, j, nnodes)] != PRESENT) || (a[CMC(j, i, nnodes)] != PRESENT))
        continue;

      if (c_directed_path(i, j, a, nnodes, nodes, path, scratch, 0)) {

        if (vstructures_and_cycles(a, nnodes, i, j)) {

          a[CMC(i, j, nnodes)] = FIXED;
          a[CMC(j, i, nnodes)] = FIXED;
          if (debuglevel > 0)
            Rprintf("  > fixing arc %s - %s, %s -> %s introduces new v-structures, %s -> %s creates cycles (step 3a).\n", \
            NODE(i), NODE(j), NODE(j), NODE(i), NODE(i), NODE(j));

        }/*THEN*/
        else {

          a[CMC(i, j, nnodes)] = FIXED;
          a[CMC(j, i, nnodes)] = ABSENT;

          if (debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(i), NODE(j));

        }/*THEN*/

        changed = CHANGED;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  Free1D(path);
  Free1D(scratch);

  return changed;

}/*PREVENT_CYCLES*/

short int creates_unshielded_collider(int *a, int nnodes, int j, int i) {

int k = 0;

  /* the target arc is NODE(j) -> NODE(i); so to have an unshielded collider
   * the graph must contain NODE(k) -> NODE(i) and there must be no arc
   * between NODE(j) and NODE(K), for at least one NODE(k). */
  for (k = 0; k < nnodes; k++) {

    if ((k == i) || (k == j))
      continue;

    if ((a[CMC(k, i, nnodes)] >= PRESENT) && (a[CMC(i, k, nnodes)] == ABSENT) &&
        (a[CMC(k, j, nnodes)] == ABSENT) && (a[CMC(j, k, nnodes)] == ABSENT))

      return TRUE;

  }/*FOR*/

  return FALSE;

}/*CREATES_UNSHIELDED_COLLIDER*/

static int prevent_additional_vstructures(int *a, SEXP nodes, int nnodes,
    int debuglevel) {

int i = 0, j = 0;
int *path = NULL, *scratch = NULL;
short int *has_parent = NULL, *has_neighbour = NULL;
short int changed = FALSE, unshielded_i = FALSE, unshielded_j = FALSE;
short int path_to_i = FALSE, path_to_j = FALSE;

  has_parent = Calloc1D(nnodes, sizeof(short int));
  has_neighbour = Calloc1D(nnodes, sizeof(short int));
  path = Calloc1D(nnodes, sizeof(int));
  scratch = Calloc1D(nnodes, sizeof(int));

  for (j = 0; j < nnodes; j++) {

    /* check if this node has both a parent (which implies a directed, incoming
     * arc) and a neighbour (which implies an undirected incident arc).  */
    for (i = 0; i < nnodes; i++) {

       if (a[CMC(i, j, nnodes)] != ABSENT) {

         if (a[CMC(j, i, nnodes)] == PRESENT)
           has_neighbour[j] = TRUE;
         else
           has_parent[j] = TRUE;

       }/*THEN*/

       /* stop at once if both are present. */
       if (has_parent[j] && has_neighbour[j]) break;

    }/*FOR*/

  }/*FOR*/

  for (j = 0; j < nnodes; j++) {

    /* all neighbours must be changed into children to prevent the creation
     * of new v-structures.*/
    if (!(has_parent[j] && has_neighbour[j]))
      continue;

    for (i = 0; i < nnodes; i++) {

      /* the target arc must be directed. */
      if (!((a[CMC(i, j, nnodes)] != ABSENT) && (a[CMC(j, i, nnodes)] == PRESENT)))
        continue;

      /* both nodes incident on the arc have parents, be careful not to
       * introduce v-structures centered on either node. */
      if (has_parent[i]) {

        unshielded_i = creates_unshielded_collider(a, nnodes, j, i);
        unshielded_j = creates_unshielded_collider(a, nnodes, i, j);
        path_to_i = c_directed_path(j, i, a, nnodes, nodes, path, scratch, 0);
        path_to_j = c_directed_path(i, j, a, nnodes, nodes, path, scratch, 0);

        if (unshielded_j && !unshielded_i) {

          /* if NODE(i) -> NODE(j) introduces an unshielded collider, fix
           * NODE(j) -> NODE(i) unless it creates cycles. */
          if (path_to_j) {

            if (debuglevel > 0)
              Rprintf("  @ arc %s -> %s introduces cycles, %s -> %s introduces an unshielded collider (step 3b).\n",
                NODE(i), NODE(j), NODE(j), NODE(i));

            goto failsafe;

          }/*THEN*/

          a[CMC(i, j, nnodes)] = ABSENT;
          a[CMC(j, i, nnodes)] = FIXED;
          changed = TRUE;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, %s -> %s introduces an unshielded collider (step 3b).\n",
              NODE(i), NODE(j), NODE(j), NODE(i));

        }/*THEN*/
        else if (unshielded_j && !unshielded_i) {

          /* if NODE(j) -> NODE(i) introduces an unshielded collider, fix
           * NODE(i) -> NODE(j) unless it creates cycles. */
          if (path_to_i) {

            if (debuglevel > 0)
              Rprintf("  @ arc %s -> %s introduces cycles, %s -> %s introduces an unshielded collider (step 3b).\n",
                NODE(j), NODE(i), NODE(i), NODE(j));

            goto failsafe;

          }/*THEN*/

          a[CMC(i, j, nnodes)] = FIXED;
          a[CMC(j, i, nnodes)] = ABSENT;
          changed = TRUE;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, %s -> %s introduces an unshielded collider (step 3b).\n",
              NODE(j), NODE(i), NODE(i), NODE(j));

        }/*THEN*/
        else if (unshielded_i && unshielded_j) {

failsafe:

          /* both directions introduce an unshielded collider, so fix the arc
           * as undirected. */
          a[CMC(i, j, nnodes)] = FIXED;
          a[CMC(j, i, nnodes)] = FIXED;
          changed = TRUE;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s - %s, both nodes have parents (step 3b).\n", NODE(j), NODE(i));

        }/*THEN*/
        else {

          /* if neither direction introduces an unshielded collider, do no not
             do anything. */

        }/*ELSE*/

      }/*THEN*/
      else {

        if (!creates_unshielded_collider(a, nnodes, i, j)) {

          if (debuglevel > 0)
            Rprintf("  > arc %s - %s is not part of an unshielded collider, leave it undirected.\n", NODE(j), NODE(i));

        }/*THEN*/
        else if (c_directed_path(i, j, a, nnodes, nodes, path, scratch, 0)) {

          a[CMC(i, j, nnodes)] = FIXED;
          a[CMC(j, i, nnodes)] = FIXED;
          changed = TRUE;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s - %s, one direction introduces new v-structures, the other creates cycles (step 3b).\n",
              NODE(j), NODE(i));

        }/*THEN*/
        else {

          a[CMC(i, j, nnodes)] = ABSENT;
          a[CMC(j, i, nnodes)] = FIXED;
          changed = TRUE;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, prevent v-structure (step 3b).\n", NODE(j), NODE(i));

        }/*THEN*/

      }/*ELSE*/

    }/*FOR*/

  }/*FOR*/

  Free1D(has_neighbour);
  Free1D(has_parent);
  Free1D(path);
  Free1D(scratch);

  return changed;

}/*PREVENT_ADDITIONAL_VSTRUCTURES*/

static int creates_conflicts(int *a, int nnodes, short int *scratch) {

int k = 0, l = 0;

  for (k = 0; k < nnodes; k++) {

    if (!scratch[k])
      continue;

    for (l = k + 1; l < nnodes; l++) {

      if (!scratch[l])
        continue;

      /* fix NODE(i) -> NODE(j) if there is such a pair of nodes. */
      if ((a[CMC(k, l, nnodes)] == ABSENT) && (a[CMC(l, k, nnodes)] == ABSENT))
        return TRUE;

    }/*FOR*/

  }/*FOR*/

  return FALSE;

}/*CREATES_CONFLICTS*/

static int prevent_conflicts(int *a, SEXP nodes, int nnodes, short int *scratch,
    int debuglevel) {

int i = 0, j = 0, k = 0, changed = FALSE;

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to do, this is not a valid arc. */
      if (i == j)
        continue;

      /* nothing to do, this is not an undirected arc. */
      if (!((a[CMC(i, j, nnodes)] == PRESENT) && (a[CMC(j, i, nnodes)] == PRESENT)))
        continue;

      /* reset the scratch space. */
      memset(scratch, '\0', nnodes * sizeof(short int));

      /* look for nodes NODE(K) such that NODE(i) -- NODE(K) -> NODE(j). */
      for (k = 0; k < nnodes; k++) {

        /* it must be a different node. */
        if ((k == i) || (k == j))
          continue;

        if ((a[CMC(i, k, nnodes)] == PRESENT) && (a[CMC(k, i, nnodes)] == PRESENT) &&
            (a[CMC(k, j, nnodes)] >= PRESENT) && (a[CMC(j, k, nnodes)] == ABSENT)) {

          scratch[k] = TRUE;

          if (debuglevel)
            Rprintf("    @ found chain %s - %s -> %s.\n",
              NODE(i), NODE(k), NODE(j));

        }/*THEN*/

      }/*FOR*/

      /* check if any pair of such nodes are not adjacent; if such a pair exists
       * fix NODE(i) -> NODE(j). */
      if (creates_conflicts(a, nnodes, scratch)) {

        a[CMC(i, j, nnodes)] = FIXED;
        a[CMC(j, i, nnodes)] = ABSENT;
        changed = TRUE;

        if (debuglevel)
          Rprintf("    > fixing %s -> %s because of these chains (step 3c).\n",
            NODE(i), NODE(j));

        break;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return changed;

}/*PREVENT_CONFLICTS*/

static int find_chain(int *a, int nnodes, int *elem) {

int i = elem[0], j = elem[1], k = 0, l = 0;

  /* looking for a chain NODE(i) -- NODE(k) -> NODE(l) -> NODE(j). */
  for (k = 0; k < nnodes; k++) {

    if (!((a[CMC(i, k, nnodes)] >= PRESENT) && (a[CMC(k, i, nnodes)] >= PRESENT)))
      continue;

    for (l = 0; l < nnodes; l++) {

      if (!((a[CMC(k, l, nnodes)] >= PRESENT) && (a[CMC(l, k, nnodes)] == ABSENT) &&
            (a[CMC(l, j, nnodes)] >= PRESENT) && (a[CMC(j, l, nnodes)] == ABSENT)))
        continue;

      /* check whether NODE(i) -- NODE(k) and NODE(i) -- NODE(l) but
       * NODE(j) is not adjacent to NODE(l). */
      if ((a[CMC(i, l, nnodes)] >= PRESENT) && (a[CMC(l, i, nnodes)] >= PRESENT) &&
          (a[CMC(k, j, nnodes)] == ABSENT) && (a[CMC(j, k, nnodes)] == ABSENT)) {

        /* save the nodes in the chain, to use in the debugging output. */
        elem[0] = k;
        elem[1] = l;

        return TRUE;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return FALSE;

}/*FIND_CHAIN*/

static int prevent_chains(int *a, SEXP nodes, int nnodes, int debuglevel) {

int i = 0, j = 0, changed = FALSE, chain[2] = {0, 0};

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to do, this is not a valid arc. */
      if (i == j)
        continue;

      /* nothing to do, this is not an undirected arc. */
      if (!((a[CMC(i, j, nnodes)] == PRESENT) && (a[CMC(j, i, nnodes)] == PRESENT)))
        continue;

      chain[0] = i;
      chain[1] = j;

      if (find_chain(a, nnodes, chain)) {

        a[CMC(i, j, nnodes)] = FIXED;
        a[CMC(j, i, nnodes)] = ABSENT;
        changed = TRUE;

        if (debuglevel)
          Rprintf("    @ found chain %s - %s -> %s -> %s with %s - %s and %s - %s, "
                  "fixing %s -> %s (step 3d).\n",
            NODE(i), NODE(j), NODE(chain[0]), NODE(chain[1]),
            NODE(i), NODE(chain[0]), NODE(i), NODE(chain[1]), NODE(i), NODE(j));

        break;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return changed;

}/*PREVENT_CHAINS*/

static void renormalize_amat(int *a, int nnodes) {

  for (int i = 0; i < nnodes * nnodes; i++)
    a[i] = (a[i] > 1) ? 1 : a[i];

}/*RENORMALIZE_AMAT*/

