#include "include/rcore.h"
#include "include/graph.h"
#include "include/allocations.h"
#include "include/matrix.h"

#define ABSENT     0
#define PRESENT    1
#define FIXED      2

#define UNCHANGED  0
#define CHANGED    1

static void fix_all_directed(int *a, int *nnodes);
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
static void is_a_sink(int *a, int node, int *k, int nnodes, int *nbr,
    short int *matched);
static int all_adjacent(int *a, int node, int k, int nnodes, int *nbr);

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
  collider = allocstatus(nnodes);

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

  return result;

}/*VSTRUCTURES*/

static SEXP amat2vstructs(int *a, SEXP nodes, int *nnodes, short int *collider) {

int i = 0, j = 0, k = 0, nvstructs = 0, counter = 0, row = 0;
SEXP result, dimnames;

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
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, mkStringVec(3, "X", "Z", "Y"));
  setAttrib(result, R_DimNamesSymbol, dimnames);

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

  UNPROTECT(2);

  return result;

}/*AMAT2VSTRUCTS*/

SEXP cpdag(SEXP arcs, SEXP nodes, SEXP moral, SEXP fix, SEXP debug) {

int i = 0, changed = 0, nnodes = length(nodes);
short int *collider = NULL;
int *a = NULL, debuglevel = isTRUE(debug);
SEXP amat;

  /* build the adjacency matrix and dereference it. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* allocate and initialize a status vector to identify the colliders. */
  collider = allocstatus(nnodes);

  /* STEP 1: scan the graph. */
  scan_graph(a, nodes, &nnodes, collider, debuglevel);

  if (isTRUE(fix)) {

    /* skipping STEP 2, just propagate the directions from the arcs that are
     * diredcted (this is for learning algorithms, to preserve whitelisted and
     * blacklisted arcs). */

    fix_all_directed(a, &nnodes);

  }/*THEN*/
  else {

    /* STEP 2: all the arcs not part of a v-structure are now undirected. */
    mark_vstructures(a, nodes, &nnodes, collider, debuglevel);

    /* decide whether or not to keep moral v-structures (i.e. shielded colliders). */
    if (isTRUE(moral))
      unmark_shielded(a, nodes, &nnodes, collider, debuglevel);

  }/*ELSE*/

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

  return amat;

}/*CPDAG*/

static void fix_all_directed(int *a, int *nnodes) {

int i = 0, j = 0;

  for (i = 0; i < *nnodes; i++) {

    for (j = 0; j < *nnodes; j++) {

      if ((a[CMC(i, j, *nnodes)] == PRESENT) && (a[CMC(j, i, *nnodes)] == ABSENT))
        a[CMC(i, j, *nnodes)] = FIXED;
      if ((a[CMC(i, j, *nnodes)] == ABSENT) && (a[CMC(j, i, *nnodes)] == PRESENT))
        a[CMC(j, i, *nnodes)] = FIXED;

    }/*FOR*/

  }/*FOR*/

}/*FIX_ALL_DIRECTED*/

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

  for (j = 0; j < n; j++) {

    for (i = j + 1; i < n; i++) {

      if ((a[CMC(i, j, n)] != PRESENT) || (a[CMC(j, i, n)] != PRESENT))
        continue;

      if (c_directed_path(i, j, a, n, nodes, 0)) {

        a[CMC(i, j, n)] = FIXED;
        a[CMC(j, i, n)] = ABSENT;

        if (debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(i), NODE(j));

        changed = CHANGED;

      }/*THEN*/
      else if (c_directed_path(j, i, a, n, nodes, 0)) {

        a[CMC(i, j, n)] = ABSENT;
        a[CMC(j, i, n)] = FIXED;

        if (debuglevel)
            Rprintf("  > fixing arc %s -> %s due to directed path (step 3a).\n", NODE(j), NODE(i));

        changed = CHANGED;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return changed;

}/*PREVENT_CYCLES*/

static int prevent_additional_vstructures(int *a, SEXP nodes, int *nnodes,
    short int *collider, int debuglevel) {

int i = 0, j = 0;
short int *has_parent = NULL, *has_neighbour = NULL;

  has_parent = allocstatus(*nnodes);
  has_neighbour = allocstatus(*nnodes);

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

      return CHANGED;

    }/*THEN*/

  }/*FOR*/

  return UNCHANGED;

}/*PREVENT_ADDITIONAL_VSTRUCTURES*/

static void renormalize_amat(int *a, int *nnodes) {

  for (int i = 0; i < (*nnodes) * (*nnodes); i++)
    a[i] = (a[i] > 1) ? 1 : a[i];

}/*RENORMALIZE_AMAT*/

/* construct a consistent DAG extension of a CPDAG. */
SEXP pdag_extension(SEXP arcs, SEXP nodes, SEXP debug) {

int i = 0, j = 0, k = 0, t = 0, nnodes = length(nodes);
int changed = 0, left = nnodes;
int *a = NULL, *nbr = NULL, debuglevel = isTRUE(debug);
short int *matched = NULL;
SEXP amat, result;

  /* build and dereference the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* aqllocate and initialize the neighbours and matched vectors. */
  nbr = alloc1dcont(nnodes);
  matched = allocstatus(nnodes);

  for (t = 0; t < nnodes; t++) {

    if (debuglevel > 0) {

      Rprintf("----------------------------------------------------------------\n");
      Rprintf("> performing pass %d.\n", t + 1);
      Rprintf("> candidate nodes: ");
        for (j = 0; j < nnodes; j++)
          if (matched[j] == 0)
            Rprintf("%s ", NODE(j));
      Rprintf("\n");

    }/*THEN*/

    for (i = 0; i < nnodes; i++) {

      /* if the node is already ok, skip it. */
      if (matched[i] != 0)
        continue;

      /* check whether the node is a sink. */
      is_a_sink(a, i, &k, nnodes, nbr, matched);

      /* if the node is not a sink move on. */
      if (k == -1) {

        if (debuglevel > 0)
          Rprintf("  * node %s is not a sink.\n", NODE(i));

        continue;

      }/*THEN*/
      else {

        if (debuglevel > 0)
          Rprintf("  * node %s is a sink.\n", NODE(i));

      }/*ELSE*/

      if (!all_adjacent(a, i, k, nnodes, nbr)) {

        if (debuglevel > 0)
          Rprintf("  * not all nodes linked to %s by an undirected arc are adjacent.\n", NODE(i));

        continue;

      }/*THEN*/
      else {

        if (debuglevel > 0) {

          if (k == 0)
            Rprintf("  * no node is linked to %s by an undirected arc.\n", NODE(i));
          else
            Rprintf("  * all nodes linked to %s by an undirected arc are adjacent.\n", NODE(i));

        }/*THEN*/

      }/*ELSE*/

      /* the current node meets all the conditions, direct all the arcs towards it. */
      if (k == 0) {

        if (debuglevel > 0)
          Rprintf("  @ no undirected arc to direct towards %s.\n", NODE(i));

      }/*THEN*/
      else {

        for (j = 0; j < k; j++)
          a[CMC(i, nbr[j], nnodes)] = 0;

        if (debuglevel > 0)
          Rprintf("  @ directing all incident undirected arcs towards %s.\n", NODE(i));

      }/*ELSE*/

      /* set the changed flag. */
      changed = 1;

      /* exclude the node from later iterations. */
      matched[i] = 1;
      left--;

    }/*FOR*/

    /* if nothing changed in the last iteration or there are no more candidate
     * nodes, there is nothing else to do. */
    if ((changed == 0) || (left == 0))
      break;
    else
      changed = 0;

  }/*FOR*/

  /* build the new arc set from the adjacency matrix. */
  PROTECT(result = amat2arcs(amat, nodes));

  UNPROTECT(2);

  return result;

}/*PDAG_EXTENSION*/

static void is_a_sink(int *a, int node, int *k, int nnodes, int *nbr,
    short int *matched) {

int j = 0;

  /* check whether the current node has outgoing arcs. */
  for (j = 0, *k = 0; j < nnodes; j++) {

    if (matched[j] != 0)
      continue;

    if ((a[CMC(j, node, nnodes)] == 0) && (a[CMC(node, j, nnodes)] == 1)) {

      /* this node is not a candidate, go to the next one. */
      *k = -1;

      break;

    }/*THEN*/
    else if ((a[CMC(j, node, nnodes)] == 1) && (a[CMC(node, j, nnodes)] == 1)) {

      /* get the nodes which are linked to the current one by an
       * undirected arc. */
      nbr[(*k)++] = j;

    }/*THEN*/

  }/*FOR*/

}/*IS_A_SINK*/

static int all_adjacent(int *a, int node, int k, int nnodes, int *nbr) {

int j = 0, l = 0;

  for (j = 0; j < k; j++) {

    for (l = j + 1; l < k; l++) {

      if ((a[CMC(nbr[j], nbr[l], nnodes)] == 0) &&
          (a[CMC(nbr[l], nbr[j], nnodes)] == 0)) {

        /* this node is not a candidate, go to the next one. */
        return FALSE;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  return TRUE;

}/*ALL_ADJACENT*/


