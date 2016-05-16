#include "include/rcore.h"
#include "include/graph.h"
#include "include/globals.h"
#include "include/sampling.h"
#include "include/bn.h"
#include "include/matrix.h"

static SEXP bn_base_structure(SEXP nodes, SEXP args, SEXP arcs, SEXP cached,
    double ntests, char *test, char *algo);

static SEXP ic_2nodes(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug);

static SEXP c_ide_cozman(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug);

static int ic_logic(int *amat, SEXP nodes, int *nnodes, int *arc, int *work,
    int *degree, double *max, int *in_degree, double *max_in, int *out_degree,
    double *max_out, int *cozman, int *path, int *scratch, int debuglevel);

static void print_modelstring(SEXP bn);

/* generate an empty graph. */
SEXP empty_graph(SEXP nodes, SEXP num) {

int i = 0, nnodes = length(nodes), *n = INTEGER(num);
SEXP list, res, args, arcs, cached;
SEXP elnames, base, base2;

  /* an empty list of optional arguments. */
  PROTECT(args = allocVector(VECSXP, 0));

  /* names for the cached information. */
  PROTECT(elnames = mkStringVec(4, "mb", "nbr", "parents", "children"));

  /* allocate and initialize the arc set. */
  PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
  setDimNames(arcs, R_NilValue, mkStringVec(2, "from", "to"));

  /* allocate and initialize nodes' cached information. */
  PROTECT(base2 = allocVector(STRSXP, 0));
  PROTECT(base = allocVector(VECSXP, 4));
  setAttrib(base, R_NamesSymbol, elnames);

  PROTECT(cached = allocVector(VECSXP, nnodes));
  setAttrib(cached, R_NamesSymbol, nodes);

  for (i = 0; i < 4; i++)
    SET_VECTOR_ELT(base, i, base2);
  for (i = 0; i < nnodes; i++)
    SET_VECTOR_ELT(cached, i, base);

  /* generate the "bn" structure. */
  PROTECT(res = bn_base_structure(nodes, args, arcs, cached, 0, "none", "empty"));

  /* return a list if more than one bn is generated. */
  if (*n > 1) {

    PROTECT(list = allocVector(VECSXP, *n));
    for (i = 0; i < *n; i++)
      SET_VECTOR_ELT(list, i, res);

    UNPROTECT(8);
    return list;

  }/*THEN*/
  else {

    UNPROTECT(7);
    return res;

  }/*ELSE*/

}/*EMPTY_GRAPH*/

/* generate a graph with given node ordering and arc probability. */
SEXP ordered_graph(SEXP nodes, SEXP num, SEXP prob) {

int i = 0, j = 0, k = 0, nnodes = length(nodes), *a = NULL, *n = INTEGER(num);
double *p = REAL(prob);
SEXP list, res, args, amat, arcs, cached, null, temp;

  /* the list of optional arguments. */
  PROTECT(args = allocVector(VECSXP, 1));
  setAttrib(args, R_NamesSymbol, mkString("prob"));
  SET_VECTOR_ELT(args, 0, prob);

  /* allocate and initialize the adjacency matrix. */
  PROTECT(amat = allocMatrix(INTSXP, nnodes, nnodes));
  a = INTEGER(amat);
  memset(a, '\0', nnodes * nnodes * sizeof(int));

  GetRNGstate();

#define ORDERED_AMAT(prob) \
      for (i = 0; i < nnodes; i++) \
        for (j = i + 1; j < nnodes; j++) \
          if (unif_rand() < prob) \
            a[CMC(i, j, nnodes)] = 1; \
          else \
            a[CMC(i, j, nnodes)] = 0; \

  /* return a list if more than one bn is generated. */
  if (*n > 1) {

    PROTECT(list = allocVector(VECSXP, *n));
    PROTECT(null = allocVector(NILSXP, 1));

    /* generate the "bn" structure, with dummy NULLs for the "arcs" and
     * "nodes" elements (which will be initialized later on). */
    PROTECT(res = bn_base_structure(nodes, args, null, null, 0, "none", "ordered"));

    for (k = 0; k < *n; k++) {

      /* sample each arc in the upper-triangular portion of the adjacency matrix
       * (so that node ordering is conserved) with the specified probability. */
      ORDERED_AMAT(*p);

      /* generate the arc set and the cached information form the adjacency
       * matrix. */
      PROTECT(arcs = amat2arcs(amat, nodes));
      PROTECT(cached = cache_structure(nodes, amat, FALSESEXP));
      SET_VECTOR_ELT(res, 1, cached);
      SET_VECTOR_ELT(res, 2, arcs);

      /* save the structure in the list. */
      PROTECT(temp = duplicate(res));
      SET_VECTOR_ELT(list, k, temp);

      UNPROTECT(3);

    }/*FOR*/

    PutRNGstate();

    UNPROTECT(5);
    return list;

  }/*THEN*/
  else {

    /* sample each arc in the upper-triangular portion of the adjacency matrix
     * (so that node ordering is conserved) with the specified probability. */
    ORDERED_AMAT(*p);

    /* generate the arc set and the cached information form the adjacency
     * matrix. */
    PROTECT(arcs = amat2arcs(amat, nodes));
    PROTECT(cached = cache_structure(nodes, amat, FALSESEXP));

    /* generate the "bn" structure. */
    PROTECT(res = bn_base_structure(nodes, args, arcs, cached, 0, "none", "ordered"));

    PutRNGstate();

    UNPROTECT(5);
    return res;

  }/*ELSE*/

}/*ORDERED_GRAPH*/

/* generate a connected graph with uniform probability, subject to some
 * constraints on the degree of the nodes. */
SEXP ide_cozman_graph(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug) {

SEXP graphlist;

  switch(length(nodes)) {

    case 1:

      /* there is only one graph with 1 node, and it's empty. */
      graphlist = empty_graph(nodes, num);
      break;

    case 2:
      /* Ide-Cozman has no mixing with only 2 nodes, work around with i.i.d.
       * sampling.*/
      if (isTRUE(connected)) {

        graphlist = ic_2nodes(nodes, num, burn_in, max_in_degree,
                      max_out_degree, max_degree, connected, debug);
        break;

      }/*THEN*/

    default:

      graphlist = c_ide_cozman(nodes, num, burn_in, max_in_degree,
                    max_out_degree, max_degree, connected, debug);

  }/*SWITCH*/

  return graphlist;

}/*IDE_COZMAN_GRAPH*/

/* an Ide-Cozman alternative for 2-nodes graphs. */
static SEXP ic_2nodes(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug) {

int i = 0, *n = INTEGER(num), *a = NULL;
int debuglevel = isTRUE(debug);
double u = 0;
SEXP list, resA, resB, arcsA, arcsB, cachedA, cachedB;
SEXP amatA, amatB, args;

  /* the list of optional arguments. */
  PROTECT(args = allocVector(VECSXP, 4));
  setAttrib(args, R_NamesSymbol,
    mkStringVec(4, "burn.in", "max.in.degree", "max.out.degree", "max.degree"));
  SET_VECTOR_ELT(args, 0, burn_in);
  SET_VECTOR_ELT(args, 1, max_in_degree);
  SET_VECTOR_ELT(args, 2, max_out_degree);
  SET_VECTOR_ELT(args, 3, max_degree);

  /* allocate and initialize the tow adjacency matrices. */
  PROTECT(amatA = allocMatrix(INTSXP, 2, 2));
  a = INTEGER(amatA);
  memset(a, '\0', sizeof(int) * 4);
  a[2] = 1;
  PROTECT(amatB = allocMatrix(INTSXP, 2, 2));
  a = INTEGER(amatB);
  memset(a, '\0', sizeof(int) * 4);
  a[1] = 1;
  /* generates the arc sets. */
  PROTECT(arcsA = amat2arcs(amatA, nodes));
  PROTECT(arcsB = amat2arcs(amatB, nodes));
  /* generate the cached node information. */
  PROTECT(cachedA = cache_structure(nodes, amatA, FALSESEXP));
  PROTECT(cachedB = cache_structure(nodes, amatB, FALSESEXP));
  /* generate the two "bn" structures. */
  PROTECT(resA = bn_base_structure(nodes, args, arcsA, cachedA, 0, "none", "empty"));
  PROTECT(resB = bn_base_structure(nodes, args, arcsB, cachedB, 0, "none", "empty"));

  if (debuglevel > 0)
    Rprintf("* no burn-in required.\n");

  GetRNGstate();

  /* return a list if more than one bn is generated. */
  if (*n > 1) {

    PROTECT(list = allocVector(VECSXP, *n));
    for (i = 0; i < *n; i++) {

      if (debuglevel > 0)
        Rprintf("* current model (%d):\n", i + 1);

      /* sample which graph to return. */
      u = unif_rand();

      if (u <= 0.5) {

        /* pick the graph with A -> B. */
        SET_VECTOR_ELT(list, i, resA);

        /* print the model string to allow a sane debugging experience. */
        if (debuglevel > 0)
          print_modelstring(resA);

      }/*THEN*/
      else {

        /* pick the graph with B -> A. */
        SET_VECTOR_ELT(list, i, resB);

        /* print the model string to allow a sane debugging experience. */
        if (debuglevel > 0)
          print_modelstring(resB);

      }/*ELSE*/

    }/*FOR*/

    PutRNGstate();

    UNPROTECT(10);
    return list;

  }/*THEN*/
  else {

    if (debuglevel > 0)
      Rprintf("* current model (1):\n");

    /* sample which graph to return. */
    u = unif_rand();

    PutRNGstate();

    UNPROTECT(9);

    if (u <= 0.5) {

      /* print the model string to allow a sane debugging experience. */
      if (debuglevel > 0)
        print_modelstring(resA);

      /* return the graph with A -> B. */
      return resA;

    }/*THEN*/
    else {

      /* print the model string to allow a sane debugging experience. */
      if (debuglevel > 0)
        print_modelstring(resB);

      /* return the graph with B -> A. */
      return resB;

    }/*ELSE*/

  }/*ELSE*/

}/*IC_2NODES*/

static SEXP c_ide_cozman(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug) {

int i = 0, k = 0, nnodes = length(nodes), *n = INTEGER(num);
int changed = 0, *work = NULL, *arc = NULL, *a = NULL, *burn = INTEGER(burn_in);
int *degree = NULL, *in_degree = NULL, *out_degree = NULL;
int *path = NULL, *scratch = NULL;
int debuglevel = isTRUE(debug), *cozman = LOGICAL(connected);
double *max_in = REAL(max_in_degree), *max_out = REAL(max_out_degree),
  *max = REAL(max_degree);
SEXP list, res, args, amat, arcs, cached, null, temp;
char *label = (*cozman > 0) ? "ic-dag" : "melancon";

  /* the list of optional arguments. */
  PROTECT(args = allocVector(VECSXP, 4));
  setAttrib(args, R_NamesSymbol,
    mkStringVec(4, "burn.in", "max.in.degree", "max.out.degree", "max.degree"));
  SET_VECTOR_ELT(args, 0, burn_in);
  SET_VECTOR_ELT(args, 1, max_in_degree);
  SET_VECTOR_ELT(args, 2, max_out_degree);
  SET_VECTOR_ELT(args, 3, max_degree);

  /* allocate and initialize the adjacency matrix. */
  PROTECT(amat = allocMatrix(INTSXP, nnodes, nnodes));
  a = INTEGER(amat);
  memset(a, '\0', nnodes * nnodes * sizeof(int));

  /* initialize a simple ordered tree with n nodes, where all nodes
   * have just one parent, except the first one that does not have
   * any parent. */
  for (i = 1; i < nnodes; i++)
    a[CMC(i - 1, i, nnodes)] = 1;

  /* allocate the arrays needed by SampleNoReplace. */
  arc = Calloc1D(2, sizeof(int));
  work = Calloc1D(nnodes, sizeof(int));

  /* allocate buffers for c_has_path(). */
  path = Calloc1D(nnodes, sizeof(int));
  scratch = Calloc1D(nnodes, sizeof(int));

  /* allocate and initialize the degree arrays. */
  degree = Calloc1D(nnodes, sizeof(int));
  in_degree = Calloc1D(nnodes, sizeof(int));
  out_degree = Calloc1D(nnodes, sizeof(int));

  for (i = 0; i < nnodes; i++) {

    in_degree[i] = out_degree[i] = 1;
    degree[i] = 2;

  }/*FOR*/
  in_degree[0] = out_degree[nnodes - 1] = 0;
  degree[0] = degree[nnodes - 1] = 1;

  GetRNGstate();

  /* wait for the markov chain monte carlo simulation to reach stationarity. */
  for (k = 0; k < *burn; k++) {

    if (debuglevel > 0)
      Rprintf("* current model (%d):\n", k + 1);

    changed = ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree, max_in,
                out_degree, max_out, cozman, path, scratch, debuglevel);

    /* print the model string to allow a sane debugging experience; note that this
     * has a huge impact on performance, so use it with care. */
    if ((debuglevel > 0) && (changed)) {

      PROTECT(null = allocVector(NILSXP, 1));
      PROTECT(res = bn_base_structure(nodes, args, null, null, 0, "none", label));
      PROTECT(arcs = amat2arcs(amat, nodes));
      PROTECT(cached = cache_structure(nodes, amat, FALSESEXP));
      SET_VECTOR_ELT(res, 1, cached);
      SET_VECTOR_ELT(res, 2, arcs);
      print_modelstring(res);
      UNPROTECT(4);

    }/*THEN*/

  }/*FOR*/

#define UPDATE_NODE_CACHE(cur) \
          if (debuglevel > 0) \
            Rprintf("  > updating cached information about node %s.\n", NODE(cur)); \
          memset(work, '\0', nnodes * sizeof(int)); \
          PROTECT(temp = cache_node_structure(cur, nodes, a, nnodes, work, FALSE)); \
          SET_VECTOR_ELT(cached, cur, temp); \
          UNPROTECT(1);

  /* return a list if more than one bn is generated. */
  if (*n > 1) {

    if (debuglevel > 0)
      Rprintf("* end of the burn-in iterations.\n");

    PROTECT(list = allocVector(VECSXP, *n));
    PROTECT(null = allocVector(NILSXP, 1));

    /* generate the "bn" structure, with dummy NULLs for the "arcs" and
     * "nodes" elements (which will be initialized later on). */
    PROTECT(res = bn_base_structure(nodes, args, null, null, 0, "none", label));

    for (k = 0; k < *n; k++) {

      if (debuglevel > 0)
        Rprintf("* current model (%d):\n", *burn + k + 1);

      changed = ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree,
                  max_in, out_degree, max_out, cozman, path, scratch, debuglevel);

      if (changed || (k == 0)) {

        /* generate the arc set and the cached information from the adjacency
         * matrix. */
        if (k > 0) {

          /* if a complete "bn" object is available, we can retrieve the cached
           * information about the nodes from the structure stored in the last
           * iteration and update only the elements that really need it. */
          temp = VECTOR_ELT(VECTOR_ELT(list, k - 1), 1);
          PROTECT(cached = duplicate(temp));

          /* update the first sampled nodes; both of them gain/lose either
           * a parent or a child.  */
          UPDATE_NODE_CACHE(arc[0] - 1);
          UPDATE_NODE_CACHE(arc[1] - 1);

          /* all the parents of the second sampled node gain/lose a node in
           * the markov blanket (the first sampled node, which shares a child
           * with all of them). */
          for (i = 0; i < nnodes; i++) {

            if ((i != arc[0] - 1) && (a[CMC(i, arc[1] - 1, nnodes)] == 1)) {

              UPDATE_NODE_CACHE(i);

            }/*THEN*/

          }/*FOR*/

        }/*THEN*/
        else {

          PROTECT(cached = cache_structure(nodes, amat, FALSESEXP));

        }/*ELSE*/

        PROTECT(arcs = amat2arcs(amat, nodes));
        SET_VECTOR_ELT(res, 1, cached);
        SET_VECTOR_ELT(res, 2, arcs);

        /* print the model string to allow a sane debugging experience. */
        if (debuglevel > 0)
          print_modelstring(res);

        /* save the structure in the list. */
        PROTECT(temp = duplicate(res));
        SET_VECTOR_ELT(list, k, temp);

        UNPROTECT(3);

      }/*THEN*/
      else {

        /* the adjacency matrix is unchanged; so we can just copy the bayesian
         * network from the previous iteration in the k-th slot of the list. */
        SET_VECTOR_ELT(list, k, VECTOR_ELT(list, k - 1));

      }/*ELSE*/

    }/*FOR*/

    PutRNGstate();

    Free1D(path);
    Free1D(scratch);
    Free1D(arc);
    Free1D(work);
    Free1D(degree);
    Free1D(in_degree);
    Free1D(out_degree);

    UNPROTECT(5);
    return list;

  }/*THEN*/
  else {

    if (debuglevel > 0)
      Rprintf("* end of the burn-in.\n* current model (%d):\n", *burn + 1);

    ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree,
      max_in, out_degree, max_out, cozman, path, scratch, debuglevel);

    /* generate the arc set and the cached information form the adjacency
     * matrix. */
    PROTECT(arcs = amat2arcs(amat, nodes));
    PROTECT(cached = cache_structure(nodes, amat, FALSESEXP));

    /* generate the "bn" structure. */
    PROTECT(res = bn_base_structure(nodes, args, arcs, cached, 0, "none", label));

    /* print the model string to allow a sane debugging experience. */
    if (debuglevel > 0)
      print_modelstring(res);

    PutRNGstate();

    UNPROTECT(5);

    Free1D(path);
    Free1D(scratch);
    Free1D(arc);
    Free1D(work);
    Free1D(degree);
    Free1D(in_degree);
    Free1D(out_degree);

    return res;

  }/*ELSE*/

}/*C_IDE_COZMAN*/

/* helper function which does the dirty work. */
static SEXP bn_base_structure(SEXP nodes, SEXP args, SEXP arcs, SEXP cached,
    double ntests, char *test, char *algo) {

SEXP res, learning;

  /* allocate and initialize the "learning" element. */
  PROTECT(learning = allocVector(VECSXP, 6));
  setAttrib(learning, R_NamesSymbol,
    mkStringVec(6, "whitelist", "blacklist", "test", "ntests", "algo", "args"));
  SET_VECTOR_ELT(learning, 2, mkString(test));
  SET_VECTOR_ELT(learning, 3, ScalarReal(ntests));
  SET_VECTOR_ELT(learning, 4, mkString(algo));
  SET_VECTOR_ELT(learning, 5, args);

  /* allocate and initialize the main structure. */
  PROTECT(res = allocVector(VECSXP, 3));
  setAttrib(res, R_NamesSymbol,
    mkStringVec(3, "learning", "nodes", "arcs"));
  setAttrib(res, R_ClassSymbol, mkString("bn"));
  SET_VECTOR_ELT(res, 0, learning);
  SET_VECTOR_ELT(res, 1, cached);
  SET_VECTOR_ELT(res, 2, arcs);

  UNPROTECT(2);
  return res;

}/*BN_BASE_STRUCTURE*/

static int ic_logic(int *amat, SEXP nodes, int *nnodes, int *arc, int *work,
    int *degree, double *max, int *in_degree, double *max_in, int *out_degree,
    double *max_out, int *cozman, int *path, int *scratch, int debuglevel) {

int path_exists = 0;

  /* sample an arc (that is, two nodes different from each other). */
  SampleNoReplace(2, *nnodes, arc, work);

  if (amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] == 1) {

    /* if the arc (i, j) exists in the actual graph, delete the arc,
     * provided that the underlying graph remains connected. */
    if (debuglevel > 0)
      Rprintf("  > arc %s -> %s is present.\n",
        NODE(arc[0] - 1), NODE(arc[1] - 1));

    /* skip the check on connectedness for Melancon's algorithm, it's not
     * required.  */
    if (*cozman == 1)  {

      /* if there is a(n undirected) path in the (underlying undirected) graph
       * from arc[0] to arc[1] other than arc[0] -> arc[1], the graph is still
       * connected. */
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 0;
      path_exists = c_has_path(arc[0] - 1, arc[1] - 1, amat, *nnodes, nodes,
               TRUE, FALSE, path, scratch, FALSE);
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 1;

    }/*THEN*/
    else {

      path_exists = 1;

    }/*ELSE*/

    if (path_exists) {

      if (debuglevel > 0)
        Rprintf("  @ removing arc %s -> %s.\n",
          NODE(arc[0] - 1), NODE(arc[1] - 1));

      /* update the adjacency matrix. */
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 0;

      /* update the {in,out,}degree counters. */
      in_degree[arc[1] - 1]--;
      out_degree[arc[0] - 1]--;
      degree[arc[0] - 1]--;
      degree[arc[1] - 1]--;

      return 1;

    }/*THEN*/
    else {

      if (debuglevel > 0)
        Rprintf("  @ not removing arc %s -> %s (graph not connected).\n",
          NODE(arc[0] - 1), NODE(arc[1] - 1));

      return 0;

    }/*ELSE*/

  }/*THEN*/
  else {

    /* add the arc, provided that the underlying graph remains acyclic. */
    if (debuglevel > 0)
      Rprintf("  > arc %s -> %s is not present.\n", NODE(arc[0] - 1), NODE(arc[1] - 1));

    /* do not add the arc if this violates the constraints on the degrees of the
     * nodes it's incident on. */
    if ((degree[arc[0] - 1] >= *max) || (degree[arc[1] - 1] >= *max) ||
        (out_degree[arc[0] - 1] >= *max_out) || (in_degree[arc[1] - 1] >= *max_in)) {

      if (debuglevel > 0) {

        if (degree[arc[0] - 1] >= *max)
          Rprintf("  > node %s already has degree %d, max is %lf.\n",
            NODE(arc[0] - 1), degree[arc[0] - 1], *max);
        if (degree[arc[1] - 1] >= *max)
          Rprintf("  > node %s already has degree %d, max is %lf.\n",
            NODE(arc[1] - 1), degree[arc[1] - 1], *max);
        if (out_degree[arc[0] - 1] >= *max_out)
          Rprintf("  > node %s already has out-degree %d, max is %lf.\n",
            NODE(arc[0] - 1), out_degree[arc[0] - 1], *max_out);
        if (in_degree[arc[1]] >= *max_in)
          Rprintf("  > node %s already has in-degree %d, max is %lf.\n",
            NODE(arc[1] - 1), in_degree[arc[1] - 1], *max_in);

        Rprintf("  > not adding arc %s -> %s (constraints!).\n",
          NODE(arc[0] - 1), NODE(arc[1] - 1));

      }/*THEN*/

      return 0;

    }/*THEN*/

    /* if there is a (directed) path from arc[1] to arc[0], adding arc[0] -> arc[1]
     * would create a cycle, so do not do it. */
    path_exists = c_has_path(arc[1] - 1, arc[0] - 1, amat, *nnodes,
                    nodes, FALSE, FALSE, path, scratch, FALSE);

    if (!path_exists) {

      if (debuglevel > 0)
        Rprintf("  @ adding arc %s -> %s.\n", NODE(arc[0] - 1), NODE(arc[1] - 1));

      /* update the adjacency matrix. */
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 1;

      /* update the {in,out,}degree counters. */
      in_degree[arc[1] - 1]++;
      out_degree[arc[0] - 1]++;
      degree[arc[0] - 1]++;
      degree[arc[1] - 1]++;

      return 1;

    }/*THEN*/
    else {

      if (debuglevel > 0)
        Rprintf("  > not adding arc %s -> %s (cycles!).\n",
          NODE(arc[0] - 1), NODE(arc[1] - 1));

      return 0;

    }/*ELSE*/

  }/*ELSE*/

}/*IC_LOGIC*/

/* print the model string for debugging purposes. */
static void print_modelstring(SEXP bn) {

SEXP s, t;

  /* allocate and populate the pairlist to be valuated. */
  PROTECT(t = s = allocList(2));
  SET_TYPEOF(s, LANGSXP);
  /* first slot, the function name. */
  SETCAR(t, BN_ModelstringSymbol);
  t = CDR(t);
  /* second slot, the bayesian network (the only argument). */
  SETCAR(t,  bn);
  /* evaluate ... */
  PROTECT(t = eval(s, R_GlobalEnv));
  /* ... and print the result. */
  Rprintf("  > model string is:\n%s\n", CHAR(STRING_ELT(t, 0)));
  UNPROTECT(2);

}/*PRINT_MODELSTRING*/

