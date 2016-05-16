#include "include/rcore.h"
#include "include/graph.h"
#include "include/matrix.h"

#define ABSENT     0
#define PRESENT    1
#define FIXED      2

#define UNCHANGED  0
#define CHANGED    1

static void scan_graph(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel);
static void mark_vstructures(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel);
static void unmark_shielded(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel);
static int prevent_cycles(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel);
static int prevent_additional_vstructures(int *a, SEXP nodes,
    int *nnodes, short int *collider, int debuglevel);
static void renormalize_amat(int *a, int *nnodes);
static SEXP amat2vstructs(int *a, SEXP nodes, int *nnodes, short int *collider);

/* return the v-structures present in the graph. */
SEXP vstructures(SEXP arcs, SEXP nodes, SEXP return_arcs, SEXP moral, SEXP debug) {

int i = 0, nnodes = length(nodes);
int *a = NULL, *ret_arcs = LOGICAL(return_arcs), debuglevel = isTRUE(debug);
int all_vstructs = isTRUE(moral);
short int *collider = NULL;
SEXP amat, result;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = Calloc1D(nnodes, sizeof(short int));

  /* STEP 1: scan the graph. */
  scan_graph(a, nodes, &nnodes, collider, debuglevel);

  /* STEP 2: all the arcs not part of a v-structure are now undirected. */
  mark_vstructures(a, nodes, &nnodes, collider, debuglevel);

  /* decide whether or not to keep moral v-structures (i.e. shielded colliders). */
  if (all_vstructs == 0)
    unmark_shielded(a, nodes, &nnodes, collider, debuglevel);

  /* remove all non-fixed arcs (including compelled ones). */
  for (i = 0; i < nnodes * nnodes; i++)
    a[i] = (a[i] == FIXED) ? 1 : 0;

  if (*ret_arcs > 0)
    PROTECT(result = amat2arcs(amat, nodes));
  else
    PROTECT(result = amat2vstructs(a, nodes, &nnodes, collider));

  UNPROTECT(2);

  Free1D(collider);

  return result;

}/*VSTRUCTURES*/

static SEXP amat2vstructs(int *a, SEXP nodes, int *nnodes, short int *collider) {

int i = 0, j = 0, k = 0, nvstructs = 0, counter = 0, row = 0;
SEXP result;

  for (i = 0; i < *nnodes; i++) {

    /* this is not a collider, nothing to do. */
    if (collider[i] == 0)
      continue;

    /* reset the parents' counter. */
    counter = 0;

    /* count the parents. */
    for (j = 0; j < *nnodes; j++)
      counter += a[CMC(j, i, *nnodes)];

    /* only arcs that are part of a v-structure are still present in the
     * adjacency matrix, so computing the number of v-structures is easy. */
    nvstructs += counter * (counter - 1) / 2;

  }/*FOR*/

  PROTECT(result = allocMatrix(STRSXP, nvstructs, 3));

  /* allocate and the colnames. */
  setDimNames(result, R_NilValue, mkStringVec(3, "X", "Z", "Y"));

  for (i = 0; i < *nnodes; i++) {

    /* this is not a collider, nothing to do. */
    if (collider[i] == 0)
      continue;

    for (j = 0; j < *nnodes; j++) {

      /* this is no parent. */
      if (a[CMC(j, i, *nnodes)] == 0)
        continue;

      /* this is a parent, look for the other one. */
      for (k = j + 1; k < *nnodes; k++) {

        /* again this is no parent. */
        if (a[CMC(k, i, *nnodes)] == 0)
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
              Rprintf("* fixing arc %s -> %s from the whitelist.\n", NODE(i), NODE(j));

            a[CMC(i, j, nnodes)] = FIXED;

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
              Rprintf("* fixing arc %s -> %s from the blacklist.\n", NODE(i), NODE(j));

            a[CMC(i, j, nnodes)] = FIXED;

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
short int *collider = NULL;
int *a = NULL, debuglevel = isTRUE(debug), all_fixed = FALSE;
SEXP amat;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = Calloc1D(nnodes, sizeof(short int));

  /* STEP 1a: scan the graph. */
  scan_graph(a, nodes, &nnodes, collider, debuglevel);

  /* STEP 1b: fix arcs whose direction is determined by parametric assumptions
   * (the illegal argument works like a blacklist). */
  if (illegal != R_NilValue) {

   if (debuglevel > 0)
     Rprintf("* setting the directions determined by parametric assumptions (step 1b).\n", i);

    fix_wlbl(a, nnodes, nodes, R_NilValue, illegal, FALSE);

  }/*THEN*/

  /* STEP 1c: fix arcs according to the whitelist and the blacklist. */
  if (isTRUE(fix))
    all_fixed = fix_wlbl(a, nnodes, nodes, whitelist, blacklist, debuglevel);

  if (!all_fixed) {

    /* STEP 2: all the arcs not part of a v-structure are now undirected. */
    mark_vstructures(a, nodes, &nnodes, collider, debuglevel);

    /* decide whether or not to keep moral v-structures (i.e. shielded colliders). */
    if (isTRUE(moral))
      unmark_shielded(a, nodes, &nnodes, collider, debuglevel);

  }/*THEN*/

  /* STEP 3: orient some more edges without introducing cycles or new
   * v-structures in the graph. */
  for (i = 0; i < nnodes * nnodes; i++) {

    if (debuglevel > 0)
      Rprintf("* setting the direction of more arcs (step 3, iteration %d).\n", i);

    /* STEP 3a: prevent the creation of cycles. */
    changed |= prevent_cycles(a, nodes, &nnodes, collider, debuglevel);

    /* STEP 3b: prevent the creation of additional v-structures. */
    changed |= prevent_additional_vstructures(a, nodes, &nnodes, collider, debuglevel);

    if (!changed) {

      if (debuglevel > 0)
        Rprintf("  > the graph is unchanged, skipping to the next step.\n");

      break;

    }/*THEN*/

    /* reset changed for the next iteration. */
   changed = 0;

  }/*FOR*/

  /* STEP 4: */

  /* renormalize the adjacency matrix (i.e. set all elements either to 0 or 1). */
  renormalize_amat(a, &nnodes);

  UNPROTECT(1);

  Free1D(collider);

  return amat;

}/*CPDAG*/

static void scan_graph(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel) {

 int i = 0, j = 0, counter = 0;

  /* count the parents of each node (the non-symmetric 1s in the
   * corresponding column of the adjacency matrix). */
  if (debuglevel > 0)
    Rprintf("* scanning the graph (step 1).\n");

  for (j = 0; j < *nnodes; j++) {

    counter = 0;

    for (i = 0; i < *nnodes; i++) {

      /* loops are never present in a (partially) directed acyclic graph. */
      if (i == j) continue;

      if ((a[CMC(i, j, *nnodes)] == PRESENT) && (a[CMC(j, i, *nnodes)] == ABSENT))
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

static void mark_vstructures(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel) {

 int i = 0, j = 0;

  /* set arcs belonging to a v-structure as FIXED. */
  if (debuglevel > 0)
    Rprintf("* marking v-structures (step 2).\n");

  for (j = 0; j < *nnodes; j++) {

    /* no v-structure here, mark all arcs as undirected and skip to
     * the next node. */
    if (!collider[j]) {

      for (i = 0; i < *nnodes; i++) {

        if (a[CMC(i, j, *nnodes)] == PRESENT)
          a[CMC(j, i, *nnodes)] = PRESENT;

      }/*FOR*/

    }/*THEN*/

  }/*FOR*/

  for (j = 0; j < *nnodes; j++) {

    /* no v-structure here, this node will be examined in step 3. */
    if (!collider[j]) continue;

    for (i = 0; i < *nnodes; i++) {

      /* loops are not allowed in (partially) directed graphs. */
      if (i == j) continue;

      if (a[CMC(i, j, *nnodes)] == PRESENT) {

        if (a[CMC(j, i, *nnodes)] == ABSENT) {

          /* this is a directed arc, it's part of a v-structure; mark it as FIXED. */
          a[CMC(i, j, *nnodes)] = FIXED;

          if (debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, it's part of a v-structure.\n", NODE(i), NODE(j));

        }/*THEN*/
        else if (a[CMC(j, i, *nnodes)] == PRESENT) {

          if (collider[i]) {

            /* both directions create additional v-structures, which is forbidden;
             * fix the arc without setting a direction. */
            a[CMC(i, j, *nnodes)] = a[CMC(j, i, *nnodes)] = FIXED;

            if (debuglevel > 0)
              Rprintf("  > fixing arc %s - %s due to conflicting v-structures.\n", NODE(j), NODE(i));

            warning("conflicting v-structures, the PDAG spans more than one equivalence class.");

          }/*THEN*/

        }/*THEN*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

}/*MARK_VSTRUCTURES*/

static void unmark_shielded(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel) {

int i = 0, j = 0, k = 0, check = FALSE;

  if (debuglevel > 0)
    Rprintf("* removing moral v-structures aka shielded colliders.\n");

  for (i = 0; i < *nnodes; i++) {

    /* that's not a collider, nothing to do.*/
    if (collider[i] == 0)
      continue;

    for (j = 0; j < *nnodes; j++) {

      /* this arc is not part of a v-structure. */
      if (a[CMC(j, i, *nnodes)] != FIXED)
        continue;

      if (debuglevel > 0)
        Rprintf("  > considering arc %s -> %s.\n", NODE(j), NODE(i));

      /* this is set to TRUE is the arc is part of an unshielded collider. */
      check = FALSE;

      for (k = 0; k < *nnodes; k++) {

        /* this arc is not part of a v-structure. */
        if ((a[CMC(k, i, *nnodes)] != FIXED) || (j == k))
          continue;

        if (debuglevel > 0)
          Rprintf("    > considering v-structure %s -> %s <- %s.\n", NODE(j), NODE(i), NODE(k));

        /* parents are not connected, that's a real v-structure. */
        if ((a[CMC(j, k, *nnodes)] == ABSENT) && (a[CMC(k, j, *nnodes)] == ABSENT)) {

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

        /* this arc should not be marked as FIXED and should not be directed. */
        a[CMC(j, i, *nnodes)] = a[CMC(i, j, *nnodes)] = PRESENT;

      }/*ELSE*/

    }/*FOR*/

  }/*FOR*/

}/*UNMARK_SHIELDED*/

static int prevent_cycles(int *a, SEXP nodes, int *nnodes, short int *collider,
    int debuglevel) {

int i = 0, j = 0, changed = UNCHANGED, n = *nnodes;
int *path = NULL, *scratch = NULL;

  /* allocate buffers for c_directed_path(). */
  path = Calloc1D(n, sizeof(int));
  scratch = Calloc1D(n, sizeof(int));

  for (j = 0; j < n; j++) {

    for (i = j + 1; i < n; i++) {

      if ((a[CMC(i, j, n)] != PRESENT) || (a[CMC(j, i, n)] != PRESENT))
        continue;

      if (c_directed_path(i, j, a, n, nodes, path, scratch, 0)) {

        a[CMC(i, j, n)] = FIXED;
        a[CMC(j, i, n)] = ABSENT;

        if (debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(i), NODE(j));

        changed = CHANGED;

      }/*THEN*/
      else if (c_directed_path(j, i, a, n, nodes, path, scratch, 0)) {

        a[CMC(i, j, n)] = ABSENT;
        a[CMC(j, i, n)] = FIXED;

        if (debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(j), NODE(i));

        changed = CHANGED;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  Free1D(path);
  Free1D(scratch);

  return changed;

}/*PREVENT_CYCLES*/

static int prevent_additional_vstructures(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel) {

int i = 0, j = 0;
short int *has_parent = NULL, *has_neighbour = NULL;

  has_parent = Calloc1D(*nnodes, sizeof(short int));
  has_neighbour = Calloc1D(*nnodes, sizeof(short int));

  for (j = 0; j < *nnodes; j++) {

    /* check if this node has both a parent (which implies a directed, incoming
     * arc) and a neighbour (which implies an undirected incident arc).  */
    for (i = 0; i < *nnodes; i++) {

       if (a[CMC(i, j, *nnodes)] != ABSENT) {

         if (a[CMC(j, i, *nnodes)] == PRESENT)
           has_neighbour[j] = TRUE;
         else
           has_parent[j] = TRUE;

       }/*THEN*/

       /* stop at once if both are present. */
       if (has_parent[j] && has_neighbour[j]) break;

    }/*FOR*/

  }/*FOR*/

  for (j = 0; j < *nnodes; j++) {

    /* all neighbours must be changed into children to prevent the creation
     * of new v-structures.*/
    if (has_parent[j] && has_neighbour[j]) {

      for (i = 0; i < *nnodes; i++) {

        if ((a[CMC(i, j, *nnodes)] != ABSENT) && (a[CMC(j, i, *nnodes)] == PRESENT)) {

          if (has_parent[i]) {

            a[CMC(i, j, *nnodes)] = FIXED;
            a[CMC(j, i, *nnodes)] = FIXED;

            if (debuglevel > 0)
              Rprintf("  > fixing arc %s - %s, both nodes have parents (step 3b).\n", NODE(j), NODE(i));

          }/*THEN*/
          else {

            a[CMC(i, j, *nnodes)] = ABSENT;
            a[CMC(j, i, *nnodes)] = FIXED;

            if (debuglevel > 0)
              Rprintf("  > fixing arc %s -> %s, prevent v-structure (step 3b).\n", NODE(j), NODE(i));

          }/*ELSE*/

        }/*THEN*/

      }/*FOR*/

      Free1D(has_neighbour);
      Free1D(has_parent);

      return CHANGED;

    }/*THEN*/

  }/*FOR*/

  Free1D(has_neighbour);
  Free1D(has_parent);

  return UNCHANGED;

}/*PREVENT_ADDITIONAL_VSTRUCTURES*/

static void renormalize_amat(int *a, int *nnodes) {

  for (int i = 0; i < (*nnodes) * (*nnodes); i++)
    a[i] = (a[i] > 1) ? 1 : a[i];

}/*RENORMALIZE_AMAT*/

