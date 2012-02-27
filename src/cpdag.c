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

/* return the v-structures present in the graph. */
SEXP vstructures(SEXP bn, SEXP arcs, SEXP debug) {

int *return_arcs = LOGICAL(arcs), *debuglevel = LOGICAL(debug), *nparents = NULL;
int i = 0, j = 0, k = 0, row = 0, nnodes = 0, nvstructs = 0;
SEXP nodes, node_data, parents, result, dimnames, colnames;

  /* get the nodes' data. */
  node_data = getListElement(bn, "nodes");
  nnodes = LENGTH(node_data);
  nodes = getAttrib(node_data, R_NamesSymbol);

  /* allocate and initialize the parents' counters. */
  nparents = alloc1dcont(nnodes);

  /* count the v-structures. */
  for (i = 0; i < nnodes; i++) {

    /* get the parents of the nodes. */
    parents = getListElement(VECTOR_ELT(node_data, i), "parents");
    nparents[i] = LENGTH(parents);

    if (*debuglevel > 0)
      Rprintf("* %d parents found for %s, implying %d v-structure(s).\n",
        nparents[i], NODE(i), nparents[i] * (nparents[i] - 1) / 2);

    /* count how many v-structures are centered on this node, or how many
     * arcs they involve. */
    if (*return_arcs > 0)
      nvstructs += nparents[i] * (nparents[i] > 1);
    else
      nvstructs += nparents[i] * (nparents[i] - 1) / 2;

  }/*FOR*/

  if (*debuglevel > 0)
    Rprintf("* %d v-structures present in the network.\n", nvstructs);

  /* allocate and initialize the return value. */
  if (*return_arcs > 0) {

    PROTECT(result = allocMatrix(STRSXP, nvstructs, 2));

    /* set the column names. */
    finalize_arcs(result);

    for (i = 0; i < nnodes; i++) {

      /* at least two parents are needed for a v-structure to exist. */
      if (nparents[i] < 2)
        continue;

      parents = getListElement(VECTOR_ELT(node_data, i), "parents");

      for (j = 0; j < nparents[i]; j++) {

        SET_STRING_ELT(result, CMC(row, 0, nvstructs), STRING_ELT(parents, j));
        SET_STRING_ELT(result, CMC(row, 1, nvstructs), STRING_ELT(nodes, i));
        row++;

      }/*FOR*/

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else {

    PROTECT(result = allocMatrix(STRSXP, nvstructs, 3));

    /* allocate and the colnames. */
    PROTECT(dimnames = allocVector(VECSXP, 2));
    PROTECT(colnames = allocVector(STRSXP, 3));
    SET_STRING_ELT(colnames, 0, mkChar("X"));
    SET_STRING_ELT(colnames, 1, mkChar("Z"));
    SET_STRING_ELT(colnames, 2, mkChar("Y"));
    SET_VECTOR_ELT(dimnames, 1, colnames);
    setAttrib(result, R_DimNamesSymbol, dimnames);

    for (i = 0; i < nnodes; i++) {

      /* at least two parents are needed for a v-structure to exist. */
      if (nparents[i] < 2)
        continue;

      parents = getListElement(VECTOR_ELT(node_data, i), "parents");

      for (j = 0; j < nparents[i]; j++) {

        for (k = j + 1; k < nparents[i]; k++) {

          SET_STRING_ELT(result, CMC(row, 0, nvstructs), STRING_ELT(parents, j));
          SET_STRING_ELT(result, CMC(row, 1, nvstructs), STRING_ELT(nodes, i));
          SET_STRING_ELT(result, CMC(row, 2, nvstructs), STRING_ELT(parents, k));

          row++;

        }/*FOR*/

      }/*FOR*/

    }/*FOR*/

    UNPROTECT(3);

  }/*ELSE*/

  return result;

}/*VSTRUCTURES*/

SEXP pdag_extension(SEXP arcs, SEXP nodes, SEXP debug) {

int i = 0, j = 0, k = 0, l = 0, t = 0, nnodes = LENGTH(nodes);
int changed = 0, left = nnodes;
int *a = NULL, *nbr = NULL, *debuglevel = LOGICAL(debug);
short int *matched = NULL;
SEXP amat, result;

  /* build and dereference the adjacency matrix. */
  PROTECT(amat = arcs2amat(arcs, nodes));
  a = INTEGER(amat);

  /* aqllocate and initialize the neighbours and matched vectors. */
  nbr = alloc1dcont(nnodes);
  matched = allocstatus(nnodes);

  for (t = 0; t < nnodes; t++) {

    if (*debuglevel > 0) {

      Rprintf("----------------------------------------------------------------\n");
      Rprintf("> performing pass %d.\n", t + 1);
      Rprintf("> candidate nodes: ");
        for (j = 0; j < nnodes; j++)
          if (matched[j] == 0)
            Rprintf("%s ", NODE(j));
      Rprintf("\n");

    }/*THEN*/

    for (i = 0; i < nnodes; i++) {

next_node:

      /* if the node is already ok, skip it. */
      if (matched[i] != 0)
        continue;

      /* check whether the current node has outgoing arcs. */
      for (j = 0, k = 0; j < nnodes; j++) {

        if (matched[j] != 0)
          continue;

        if ((a[CMC(j, i, nnodes)] == 0) && (a[CMC(i, j, nnodes)] == 1)) {

          if (*debuglevel > 0)
            Rprintf("  * node %s is not a sink.\n", NODE(i));

          /* this node is not a candidate, go to the next one. */
          i++;
          goto next_node;

        }/*THEN*/
        else if ((a[CMC(j, i, nnodes)] == 1) && (a[CMC(i, j, nnodes)] == 1)) {

          /* get the nodes which are linked to the current one by an 
           * undirected arc. */
          nbr[k++] = j;

        }/*THEN*/

      }/*FOR*/

      if (*debuglevel > 0)
        Rprintf("  * node %s is a sink.\n", NODE(i));

      for (j = 0; j < k; j++) {

        for (l = j + 1; l < k; l++) {

          if ((a[CMC(nbr[j], nbr[l], nnodes)] == 0) &&
              (a[CMC(nbr[l], nbr[j], nnodes)] == 0)) {

            if (*debuglevel > 0)
              Rprintf("  * not all nodes linked to %s by an undirected arc are adjacent.\n", NODE(i));

            /* this node is not a candidate, go to the next one. */
            i++;
            goto next_node;

          }/*THEN*/

        }/*FOR*/

      }/*FOR*/

      if (*debuglevel > 0)
        Rprintf("  * all nodes linked to %s by an undirected arc are adjacent.\n", NODE(i));

      /* the current node meets all the conditions, direct all the arcs towards it. */
      for (j = 0; j < k; j++)
        a[CMC(i, nbr[j], nnodes)] = 0;

      if (*debuglevel > 0)
        Rprintf("  @ directing all incident arcs towards %s.\n", NODE(i));

      /* set the changed flag. */
      changed = 1;

      /* exclude the node from later iterations. */
      matched[i] = 1;
      left--;

    }/*FOR*/

    /* if nothing changed in the last iteration or there are no more candidate
     * nodes, there is nothing else to do. */
    if (changed == 0) {

      if (left > 0)
        error("unable to construct a consistent extension.");
      else
        break;

    }/*THEN*/
    else {

      changed = 0;

    }/*ELSE*/

  }/*FOR*/

  /* build the new arc set from the adjacency matrix. */
  PROTECT(result = amat2arcs(amat, nodes));

  UNPROTECT(2);

  return result;

}/*PDAG_EXTENSION*/
