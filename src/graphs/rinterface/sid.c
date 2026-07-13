#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../graphs.h"

/* binary search for: is b reachable from a)? */
static bool reach(sparse_amat anc, int a, int b) {

int lo = anc.rowptr[b], hi = anc.rowptr[b + 1] - 1;

  while (lo <= hi) {

    int mid = (lo + hi) / 2, v = anc.colidx[mid];

    if (v == a)
      return TRUE;
    else if (v < a)
      lo = mid + 1;
    else
      hi = mid - 1;

  }/*WHILE*/

  return FALSE;

}/*REACH*/

/* This is a sparse path matrix; it stores only the (sorted) set of each node's
 * ancestors, including itself from a per-node depth-first search. */
static sparse_amat ancestor_closure(sparse_amat parents, int nnodes) {

int npairs = 0, cap = 0, top = 0;
int *from = NULL, *to = NULL, *stack = NULL, *seen = NULL;
sparse_amat anc = { 0 };

  stack = Calloc1D(nnodes, sizeof(int));
  seen = Calloc1D(nnodes, sizeof(int));
  for (int k = 0; k < nnodes; k++)
    seen[k] = -1;

  /* a starting guess for the number of reachable pairs, grown as needed. */
  cap = ((nnodes > parents.narcs) ? nnodes : parents.narcs) + nnodes;
  from = Calloc1D(cap, sizeof(int));
  to = Calloc1D(cap, sizeof(int));

  for (int b = 0; b < nnodes; b++) {

    top = 0;
    seen[b] = b;
    stack[top++] = b;

    while (top > 0) {

      int u = stack[--top];

      /* u is an ancestor of b. */
      if (npairs == cap) {

        cap *= 2;
        from = Realloc1D(from, cap, sizeof(int));
        to = Realloc1D(to, cap, sizeof(int));

      }/*THEN*/
      from[npairs] = b;
      to[npairs] = u;
      npairs++;

      for (int e = parents.rowptr[u]; e < parents.rowptr[u + 1]; e++) {

        int w = parents.colidx[e];

        if (seen[w] != b) { seen[w] = b; stack[top++] = w; }

      }/*FOR*/

    }/*WHILE*/

  }/*FOR*/

  /* fill out row b with its ancestors and sort it. */
  anc = new_sparse_amat(from, to, npairs, nnodes);

  Free1D(stack);
  Free1D(seen);
  Free1D(from);
  Free1D(to);

  return anc;

}/*ANCESTOR_CLOSURE*/

/* Inlinear-time: is j d-connected to i given the set Z (inZ) in the PROPER
 * back-door true graph with the arcs i -> c removed for every child c of i
 * that lies on a proper causal path to j (c == j or c reaches j)? */
static bool bayes_ball_dconnected(sparse_amat ch, sparse_amat pa, int p, int i,
    int j, bool *inZ, sparse_amat anc, int *top, int *bottom, int *queue,
    int epoch) {

int head = 0, tail = 0;

  /* top[X] / bottom[X] hold the epoch in which X was last queued moving up /
   * down; "X visited this call" is top[X] == epoch (no per-call O(p) reset, the
   * caller just passes a fresh epoch each time). */

  /* consume the source so it is never re-expanded as an internal node. */
  top[i] = bottom[i] = epoch;

  /* back-door: the ball goes up into the parents of i. */
  for (int e = pa.rowptr[i]; e < pa.rowptr[i + 1]; e++) {

    int u = pa.colidx[e];

    if (top[u] != epoch) {

      top[u] = epoch;
      queue[tail++] = 2 * u;

    }/*THEN*/

  }/*FOR*/

  /* forward pass: the ball goes down into the children of i, except those
   * whose arc i -> c is removed in the proper back-door graph. */
  for (int e = ch.rowptr[i]; e < ch.rowptr[i + 1]; e++) {

    int c = ch.colidx[e];

    if ((c == j) || reach(anc, c, j))
      continue;

    if (bottom[c] != epoch) {

      bottom[c] = epoch;
      queue[tail++] = 2 * c + 1;

    }/*THEN*/

  }/*FOR*/

  while (head < tail) {

    int item = queue[head++], X = item / 2, dir = item % 2;

    if (dir == 0) {

      /* visited from a child (ball moving up): an unobserved node passes the
       * ball to its parents (up) and bounces it to its children (down). */
      if (!inZ[X]) {

        for (int e = pa.rowptr[X]; e < pa.rowptr[X + 1]; e++) {

          int u = pa.colidx[e];

          if (top[u] != epoch) { top[u] = epoch; queue[tail++] = 2 * u; }

        }/*FOR*/

        for (int e = ch.rowptr[X]; e < ch.rowptr[X + 1]; e++) {

          int c = ch.colidx[e];

          if (bottom[c] != epoch) { bottom[c] = epoch; queue[tail++] = 2 * c + 1; }

        }/*FOR*/

      }/*THEN*/

    }/*THEN*/
    else {

      /* visited from a parent (ball moving down): an unobserved node passes the
       * ball to its children, an observed collider bounces it to its parents. */
      if (!inZ[X]) {

        for (int e = ch.rowptr[X]; e < ch.rowptr[X + 1]; e++) {

          int c = ch.colidx[e];

          if (bottom[c] != epoch) { bottom[c] = epoch; queue[tail++] = 2 * c + 1; }

        }/*FOR*/

      }/*THEN*/
      else {

        for (int e = pa.rowptr[X]; e < pa.rowptr[X + 1]; e++) {

          int u = pa.colidx[e];

          if (top[u] != epoch) { top[u] = epoch; queue[tail++] = 2 * u; }

        }/*FOR*/

      }/*ELSE*/

    }/*ELSE*/

  }/*WHILE*/

  return (top[j] == epoch) || (bottom[j] == epoch);

}/*BAYES_BALL_DCONNECTED*/

/* structural interventional distance. */
SEXP sid(SEXP learned, SEXP true_graph, SEXP debug) {

int nnodes = 0, count = 0, epoch = 0;
int *queue = NULL, *top = NULL, *bottom = NULL;
bool *inZ = NULL, *zdesc = NULL;
bool debugging = isTRUE(debug);
char **labels = NULL;
sparse_amat children = { 0 }, parents = { 0 }, learned_parents = { 0 };
sparse_amat ancestors = { 0 };
SEXP result;

  /* set up the children and parents adjacency matrices of the true graph, and
   * the parents adjacency of the learned graph. */
  children = sparse_amat_from_SEXP(true_graph, FALSE);
  parents = sparse_amat_from_SEXP(true_graph, TRUE);
  learned_parents = sparse_amat_from_SEXP(learned, TRUE);
  nnodes = children.dim;
  labels = children.labels;
  /* set up the path matrix of the true graph. */
  ancestors = ancestor_closure(parents, nnodes);

  inZ = Calloc1D(nnodes, sizeof(bool));
  zdesc = Calloc1D(nnodes, sizeof(bool));
  top = Calloc1D(nnodes, sizeof(int));
  bottom = Calloc1D(nnodes, sizeof(int));
  queue = Calloc1D(2 * nnodes, sizeof(int));

  /* iterate over the intervention targets... */
  for (int i = 0; i < nnodes; i++) {

    /* the adjustment set is the parents of the target in the learned graph. */
    for (int k = 0; k < nnodes; k++)
      inZ[k] = FALSE;
    for (int e = learned_parents.rowptr[i]; e < learned_parents.rowptr[i + 1]; e++)
      inZ[learned_parents.colidx[e]] = TRUE;

    /* mark the nodes that have a descendant in Z; precomputed once per target so
     * the forbidden-node test below is cheaper per pair. A node has a descendant
     * in Z iff it is an ancestor of some z in Z, so just union the ancestors of
     * the nodes in Z. */
    for (int k = 0; k < nnodes; k++)
      zdesc[k] = FALSE;
    for (int e = learned_parents.rowptr[i]; e < learned_parents.rowptr[i + 1]; e++) {

      int z = learned_parents.colidx[e];

      for (int a = ancestors.rowptr[z]; a < ancestors.rowptr[z + 1]; a++)
        zdesc[ancestors.colidx[a]] = TRUE;

    }/*FOR*/

    if (debugging) {

      Rprintf("----------------------------------------------------------------\n");
      Rprintf("* checking the causal effects of intervening on node %s.\n",
        labels[i]);
      Rprintf("  > adjustment set Z, the parents of %s in the learned graph:",
        labels[i]);
      if (learned_parents.rowptr[i] == learned_parents.rowptr[i + 1])
        Rprintf(" (empty).\n");
      else {

        for (int e = learned_parents.rowptr[i]; e < learned_parents.rowptr[i + 1]; e++)
          Rprintf(" %s", labels[learned_parents.colidx[e]]);
        Rprintf(".\n");

      }/*ELSE*/

    }/*THEN*/

    /* ... and over the other nodes. */
    for (int j = 0; j < nnodes; j++) {

      bool wrong = FALSE;

      if (j == i)
        continue;

      if (debugging)
        Rprintf("  > effect of intervening on %s on the distribution of %s:\n",
          labels[i], labels[j]);

      if (inZ[j]) {

        /* the learned model says there is no effect; wrong if the true graph
         * has a directed path from the target to the other node. */
        wrong = reach(ancestors, i, j);

        if (debugging) {

          Rprintf("    > %s is in Z, removing arcs into %s makes the learned "
            "causal effect null.\n", labels[j], labels[i]);
          Rprintf("    > ancestor matrix: %sdirected path %s ~> %s in the true "
            "graph, true effect is %s.\n", wrong ? "" : "no ", labels[i],
            labels[j], wrong ? "non-null" : "null");

        }/*THEN*/

      }/*THEN*/
      else {

        int forbidden = -1;

        if (debugging)
          Rprintf("    > %s is not in Z, checking whether Z is a valid back-door "
            "adjustment set.\n", labels[j]);

        /* forbidden nodes: is there a node m on a proper causal path i ~> j that
         * has a descendant in Z? such an m is an ancestor of j that is also a
         * descendant of i, so iterate the ancestors of j directly. */
        for (int e = ancestors.rowptr[j]; e < ancestors.rowptr[j + 1]; e++) {

          int m = ancestors.colidx[e];

          if ((m != i) && zdesc[m] && reach(ancestors, i, m)) {

            forbidden = m;
            break;

          }/*THEN*/

        }/*FOR*/

        if (forbidden >= 0) {

          /* Z adjusts for a node lying on (or below) a proper causal path,
           * blocking part of the causal effect: it is not a valid set. */
          wrong = TRUE;

          if (debugging)
            Rprintf("    > forbidden-node query: %s lies on a proper causal path "
              "%s ~> %s and has a descendant in Z, Z is not a valid adjustment "
              "set.\n", labels[forbidden], labels[i], labels[j]);

        }/*THEN*/
        else {

          /* otherwise Z is invalid iff i and j are d-connected given Z in the
           * proper back-door graph. */
          wrong = bayes_ball_dconnected(children, parents, nnodes, i, j, inZ,
                    ancestors, top, bottom, queue, ++epoch);

          if (debugging) {

            Rprintf("    > forbidden-node query: no node on a causal path has a "
              "descendant in Z.\n");
            Rprintf("    > back-door query: %s and %s are %s given Z.\n",
              labels[i], labels[j], wrong ? "d-connected" : "d-separated");

          }/*THEN*/

        }/*ELSE*/

      }/*ELSE*/

      if (debugging)
        Rprintf("    @ interventional distributions are %s.\n",
          wrong ? "different" : "the same");

      if (wrong)
        count++;

    }/*FOR*/

  }/*FOR*/

  FreeSparseAMAT(children);
  FreeSparseAMAT(parents);
  FreeSparseAMAT(learned_parents);
  FreeSparseAMAT(ancestors);
  Free1D(inZ);
  Free1D(zdesc);
  Free1D(top);
  Free1D(bottom);
  Free1D(queue);

  PROTECT(result = ScalarInteger(count));
  UNPROTECT(1);

  return result;

}/*SID*/
