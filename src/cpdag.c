#include "common.h"

#define ABSENT     0
#define PRESENT    1
#define FIXED      2

#define UNCHANGED  0
#define CHANGED    1

static void scan_graph(int *a, SEXP nodes, int *nnodes, 
    short int *collider, int *debuglevel);
static void mark_vstructures(int *a, SEXP nodes, int *nnodes, 
    short int *collider, int *debuglevel);
static int prevent_cycles(int *a, SEXP nodes, int *nnodes, 
    short int *collider, int *debuglevel);
static int prevent_additional_vstructures(int *a, SEXP nodes, 
    int *nnodes, short int *collider, int *debuglevel);
static void renormalize_amat(int *a, int *nnodes);

SEXP cpdag(SEXP arcs, SEXP nodes, SEXP debug) {

int i = 0, changed = 0, nnodes = LENGTH(nodes);
short int *collider = NULL;
int *a = NULL, *debuglevel = LOGICAL(debug);
SEXP amat;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = allocstatus(nnodes);

  /* STEP 1: scan the graph. */
  scan_graph(a, nodes, &nnodes, collider, debuglevel);

  /* STEP 2: all the arcs not part of a v-structure are now undirected. */
  mark_vstructures(a, nodes, &nnodes, collider, debuglevel);

  /* STEP 3: orient some more edges without introducing cycles or new
   * v-structures in the graph. */
  for (i = 0; i < nnodes * nnodes; i++) {

    if (*debuglevel > 0)
      Rprintf("* setting the direction of more arcs (step 3, iteration %d).\n", i);

    /* STEP 3a: prevent the creation of cycles. */
    changed |= prevent_cycles(a, nodes, &nnodes, collider, debuglevel);

    /* STEP 3b: prevent the creation of additional v-structures. */
    changed |= prevent_additional_vstructures(a, nodes, &nnodes, collider, debuglevel);

    if (!changed) {

      if (*debuglevel > 0)
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

  return amat;

}/*CPDAG*/

static void scan_graph(int *a, SEXP nodes, int *nnodes, 
    short int *collider, int *debuglevel) {

 int i = 0, j = 0, counter = 0;

  /* count the parents of each node (the non-symmetric 1s in the 
   * corresponding column of the adjacency matrix). */
  if (*debuglevel > 0)
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

        if (*debuglevel > 0)
          Rprintf("  > node %s is the center of a v-structure.\n", NODE(j));

        break;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

}/*SCAN_GRAPH*/

static void mark_vstructures(int *a, SEXP nodes, int *nnodes, 
    short int *collider, int *debuglevel) {

 int i = 0, j = 0;

  /* set arcs belonging to a v-structure as FIXED. */
  if (*debuglevel > 0)
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

          /* this is a directed arc, it's part of a v-structure; mark it as fixed. */
          a[CMC(i, j, *nnodes)] = FIXED;

          if (*debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, it's part of a v-structure.\n", NODE(i), NODE(j));

        }/*THEN*/
        else if (a[CMC(j, i, *nnodes)] == PRESENT) {

          if (collider[i]) {

            /* both directions create additional v-structures, which is forbidden; 
             * fix the arc without setting a direction. */
            a[CMC(i, j, *nnodes)] = a[CMC(j, i, *nnodes)] = FIXED;

            if (*debuglevel > 0)
              Rprintf("  > fixing arc %s - %s due to conflicting v-structures.\n", NODE(j), NODE(i));

            warning("conflicting v-structures, the PDAG spans more than one equivalence class.");

          }/*THEN*/
          else {

            /* this arc is not part of a v-structure; set its direction away from 
               the current node by exclusion.  */
            a[CMC(i, j, *nnodes)] = ABSENT;
            a[CMC(j, i, *nnodes)] = FIXED;

            if (*debuglevel > 0)
              Rprintf("  > fixing arc %s -> %s, it's not part of a v-structure.\n", NODE(j), NODE(i));

          }/*ELSE*/

        }/*THEN*/

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

}/*MARK_VSTRUCTURES*/

static int prevent_cycles(int *a, SEXP nodes, int *nnodes, short int *collider, 
    int *debuglevel) {

int i = 0, j = 0, changed = UNCHANGED, n = *nnodes;

  for (j = 0; j < n; j++) {

    for (i = j + 1; i < n; i++) {

      if ((a[CMC(i, j, n)] != PRESENT) || (a[CMC(j, i, n)] != PRESENT))
        continue;

      if (c_directed_path(i, j, a, n, nodes, 0)) {

        a[CMC(i, j, n)] = FIXED;
        a[CMC(j, i, n)] = ABSENT;

        if (*debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(j), NODE(i));

        changed = CHANGED;

      }/*THEN*/
      else if (c_directed_path(j, i, a, n, nodes, 0)) {

        a[CMC(i, j, n)] = ABSENT;
        a[CMC(j, i, n)] = FIXED;

        if (*debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(i), NODE(j));

        changed = CHANGED;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return changed;

}/*PREVENT_CYCLES*/

static int prevent_additional_vstructures(int *a, SEXP nodes, int *nnodes,
    short int *collider, int *debuglevel) {

int i = 0, j = 0;
short int has_parent = 0, has_neighbour = 0;

  for (j = 0; j < *nnodes; j++) {

    /* reset has_parent and has_neighbour. */
    has_parent = has_neighbour = 0;

    /* nodes in the middle of a v-structure were already dealt with in step 2. */
    if (collider[j] != 0) continue;

    /* check if this node has both a parent (which implies a directed, incoming
     * arc) and a neighbour (which implies an undirected incident arc).  */
    for (i = 0; i < *nnodes; i++) {

       if (a[CMC(i, j, *nnodes)] != ABSENT) {

         if (a[CMC(j, i, *nnodes)] == PRESENT)
           has_neighbour = 1;
         else
           has_parent = 1;

       }/*THEN*/

       /* stop at once if both are present. */
       if (has_parent && has_neighbour) break;

    }/*FOR*/

    /* all neighbours must be changed into children to prevent the creation
     * of new v-structures.*/
    if (has_parent && has_neighbour) {

      for (i = 0; i < *nnodes; i++) {

        if ((a[CMC(i, j, *nnodes)] != ABSENT) && (a[CMC(j, i, *nnodes)] == PRESENT)) {

          a[CMC(i, j, *nnodes)] = ABSENT;
          a[CMC(j, i, *nnodes)] = FIXED;

          if (*debuglevel > 0)
            Rprintf("  > fixing arc %s -> %s, prevent v-structure (step 3b).\n", NODE(j), NODE(i));

        }/*THEN*/

      }/*FOR*/

      return CHANGED;

    }/*THEN*/

  }/*FOR*/

  return UNCHANGED;

}/*PREVENT_ADDITIONAL_VSTRUCTURES*/

static void renormalize_amat(int *a, int *nnodes) {

  for (int i = 0; i < (*nnodes) * (*nnodes); i++)
    a[i] = (a[i] > 1) ? 1 : a[i];

}/*RENORMALIZE_AMAT*/
