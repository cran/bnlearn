#include "common.h"

static SEXP bn_base_structure(SEXP nodes, SEXP args, SEXP arcs, SEXP cached,
    double _ntests, char *_test, char *_algo);

static int ic_logic(int *amat, SEXP nodes, int *nnodes, int *arc, int *work,
    int *degree, double *max, int *in_degree, double *max_in, int *out_degree,
    double *max_out, int *cozman, int *debuglevel);

static void print_modelstring(SEXP bn);

/* generate an empty graph. */
SEXP empty_graph(SEXP nodes, SEXP num) {

int i = 0, nnodes = LENGTH(nodes), *n = INTEGER(num);
SEXP list, res, args, arcs, cached;
SEXP dimnames, colnames, elnames, base, base2;

  /* an empty list of optional arguments. */
  PROTECT(args = allocVector(VECSXP, 0));

  /* names for the arc set columns. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* names for the cached information. */
  PROTECT(elnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(elnames, 0, mkChar("mb"));
  SET_STRING_ELT(elnames, 1, mkChar("nbr"));
  SET_STRING_ELT(elnames, 2, mkChar("parents"));
  SET_STRING_ELT(elnames, 3, mkChar("children"));

  /* allocate and initialize the arc set. */
  PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
  setAttrib(arcs, R_DimNamesSymbol, dimnames);

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

    UNPROTECT(10);
    return list;

  }/*THEN*/
  else {

    UNPROTECT(9);
    return res;

  }/*ELSE*/

}/*EMPTY_GRAPH*/

/* generate a graph with given node ordering and arc probability. */
SEXP ordered_graph(SEXP nodes, SEXP num, SEXP prob) {

int i = 0, j = 0, k = 0, nnodes = LENGTH(nodes), *a = NULL, *n = INTEGER(num);
double *p = REAL(prob);
SEXP list, res, args, argnames, amat, arcs, cached, debug2, null, temp;

  /* a fake debug argument (set to FALSE) for cache_structure(). */
  PROTECT(debug2 = allocVector(LGLSXP, 1));
  LOGICAL(debug2)[0] = FALSE;

  /* the list of optional arguments. */
  PROTECT(argnames = allocVector(STRSXP, 1));
  SET_STRING_ELT(argnames, 0, mkChar("prob"));

  PROTECT(args = allocVector(VECSXP, 1));
  setAttrib(args, R_NamesSymbol, argnames);
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
      PROTECT(cached = cache_structure(nodes, amat, debug2));
      SET_VECTOR_ELT(res, 1, cached);
      SET_VECTOR_ELT(res, 2, arcs);

      /* save the structure in the list. */
      PROTECT(temp = duplicate(res));
      SET_VECTOR_ELT(list, k, temp);

      UNPROTECT(3);

    }/*FOR*/

    PutRNGstate();

    UNPROTECT(7);
    return list;

  }/*THEN*/
  else {

    /* sample each arc in the upper-triangular portion of the adjacency matrix
     * (so that node ordering is conserved) with the specified probability. */
    ORDERED_AMAT(*p);

    /* generate the arc set and the cached information form the adjacency
     * matrix. */
    PROTECT(arcs = amat2arcs(amat, nodes));
    PROTECT(cached = cache_structure(nodes, amat, debug2));

    /* generate the "bn" structure. */
    PROTECT(res = bn_base_structure(nodes, args, arcs, cached, 0, "none", "ordered"));

    PutRNGstate();

    UNPROTECT(7);
    return res;

  }/*ELSE*/

}/*ORDERED_GRAPH*/

/* generate a connected graph with uniform probability, subject to some
 * constraints on the degree of the nodes. */
SEXP ide_cozman_graph(SEXP nodes, SEXP num, SEXP burn_in, SEXP max_in_degree,
    SEXP max_out_degree, SEXP max_degree, SEXP connected, SEXP debug) {

int i = 0, k = 0, nnodes = LENGTH(nodes), *n = INTEGER(num);
int changed = 0, *work = NULL, *arc = NULL, *a = NULL, *burn = INTEGER(burn_in);
int *degree = NULL, *in_degree = NULL, *out_degree = NULL;
int *debuglevel = LOGICAL(debug), *cozman = LOGICAL(connected);
double *max_in = REAL(max_in_degree), *max_out = REAL(max_out_degree),
  *max = REAL(max_degree);
SEXP list, res, args, argnames, amat, arcs, cached, debug2, null, temp;
char *label = (*cozman > 0) ? "ic-dag" : "melancon";

  /* a fake debug argument (set to FALSE) for cache_structure(). */
  PROTECT(debug2 = allocVector(LGLSXP, 1));
  LOGICAL(debug2)[0] = FALSE;

  /* the list of optional arguments. */
  PROTECT(argnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(argnames, 0, mkChar("burn.in"));
  SET_STRING_ELT(argnames, 1, mkChar("max.in.degree"));
  SET_STRING_ELT(argnames, 2, mkChar("max.out.degree"));
  SET_STRING_ELT(argnames, 3, mkChar("max.degree"));

  PROTECT(args = allocVector(VECSXP, 4));
  setAttrib(args, R_NamesSymbol, argnames);
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
  arc = alloc1dcont(2);
  work = alloc1dcont(nnodes);

  /* allocate and initialize the degree arrays. */
  degree = alloc1dcont(nnodes);
  in_degree = alloc1dcont(nnodes);
  out_degree = alloc1dcont(nnodes);

  for (i = 0; i < nnodes; i++) {

    in_degree[i] = out_degree[i] = 1;
    degree[i] = 2;

  }/*FOR*/
  in_degree[0] = out_degree[nnodes - 1] = 0;
  degree[0] = degree[nnodes - 1] = 1;

  GetRNGstate();

  /* wait for the markov chain monte carlo simulation to reach stationarity. */
  for (k = 0; k < *burn; k++) {

    if (*debuglevel > 0)
      Rprintf("* current model (%d):\n", k + 1);

    changed = ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree, max_in,
                out_degree, max_out, cozman, debuglevel);

    /* print the model string to allow a sane debugging experience; note that this
     * has a huge impact on performance, so use it with care. */
    if ((*debuglevel > 0) && (changed)) {

      PROTECT(null = allocVector(NILSXP, 1));
      PROTECT(res = bn_base_structure(nodes, args, null, null, 0, "none", label));
      PROTECT(arcs = amat2arcs(amat, nodes));
      PROTECT(cached = cache_structure(nodes, amat, debug2));
      SET_VECTOR_ELT(res, 1, cached);
      SET_VECTOR_ELT(res, 2, arcs);
      print_modelstring(res);
      UNPROTECT(4);

    }/*THEN*/

  }/*FOR*/

#define UPDATE_NODE_CACHE(cur) \
          if (*debuglevel > 0) \
            Rprintf("  > updating cached information about node %s.\n", NODE(cur)); \
          memset(work, '\0', nnodes * sizeof(int)); \
          PROTECT(temp = c_cache_partial_structure(cur, nodes, amat, work, debug2)); \
          SET_VECTOR_ELT(cached, cur, temp); \
          UNPROTECT(1);

  /* return a list if more than one bn is generated. */
  if (*n > 1) {

    if (*debuglevel > 0)
      Rprintf("* end of the burn-in iterations.\n");

    PROTECT(list = allocVector(VECSXP, *n));
    PROTECT(null = allocVector(NILSXP, 1));

    /* generate the "bn" structure, with dummy NULLs for the "arcs" and
     * "nodes" elements (which will be initialized later on). */
    PROTECT(res = bn_base_structure(nodes, args, null, null, 0, "none", label));

    for (k = 0; k < *n; k++) {

      if (*debuglevel > 0)
        Rprintf("* current model (%d):\n", *burn + k + 1);

      changed = ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree,
                  max_in, out_degree, max_out, cozman, debuglevel);

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

          PROTECT(cached = cache_structure(nodes, amat, debug2));

        }/*ELSE*/

        PROTECT(arcs = amat2arcs(amat, nodes));
        SET_VECTOR_ELT(res, 1, cached);
        SET_VECTOR_ELT(res, 2, arcs);

        /* print the model string to allow a sane debugging experience. */
        if (*debuglevel > 0)
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

    UNPROTECT(7);
    return list;

  }/*THEN*/
  else {

    if (*debuglevel > 0)
      Rprintf("* end of the burn-in.\n* current model (%d):\n", *burn + 1);

    ic_logic(a, nodes, &nnodes, arc, work, degree, max, in_degree,
      max_in, out_degree, max_out, cozman, debuglevel);

    /* generate the arc set and the cached information form the adjacency
     * matrix. */
    PROTECT(arcs = amat2arcs(amat, nodes));
    PROTECT(cached = cache_structure(nodes, amat, debug2));

    /* generate the "bn" structure. */
    PROTECT(res = bn_base_structure(nodes, args, arcs, cached, 0, "none", label));

    /* print the model string to allow a sane debugging experience. */
    if (*debuglevel > 0)
      print_modelstring(res);

    PutRNGstate();

    UNPROTECT(7);
    return res;

  }/*ELSE*/

}/*IDE_COZMAN_GRAPH*/

/* helper function which does the dirty work. */
static SEXP bn_base_structure(SEXP nodes, SEXP args, SEXP arcs, SEXP cached,
    double _ntests, char *_test, char *_algo) {

SEXP res, learning, names, names2, class, test, ntests, algo;

  /* names of the elements of the "learning" element. */
  PROTECT(names = allocVector(STRSXP, 6));
  SET_STRING_ELT(names, 0, mkChar("whitelist"));
  SET_STRING_ELT(names, 1, mkChar("blacklist"));
  SET_STRING_ELT(names, 2, mkChar("test"));
  SET_STRING_ELT(names, 3, mkChar("ntests"));
  SET_STRING_ELT(names, 4, mkChar("algo"));
  SET_STRING_ELT(names, 5, mkChar("args"));

  /* names of the three parts of the strcture. */
  PROTECT(names2 = allocVector(STRSXP, 3));
  SET_STRING_ELT(names2, 0, mkChar("learning"));
  SET_STRING_ELT(names2, 1, mkChar("nodes"));
  SET_STRING_ELT(names2, 2, mkChar("arcs"));

  /* the class name. */
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("bn"));

  /* the number of tests/score comparisons/whatever. */
  PROTECT(ntests = allocVector(REALSXP, 1));
  NUM(ntests) = _ntests;

  /* the test used in the learning algorithm. */
  PROTECT(test = allocVector(STRSXP, 1));
  SET_STRING_ELT(test, 0, mkChar(_test));

  /* the label of the learning algorithm. */
  PROTECT(algo = allocVector(STRSXP, 1));
  SET_STRING_ELT(algo, 0, mkChar(_algo));

  /* allocate and initialize the "learning" element. */
  PROTECT(learning = allocVector(VECSXP, 6));
  setAttrib(learning, R_NamesSymbol, names);
  SET_VECTOR_ELT(learning, 2, test);
  SET_VECTOR_ELT(learning, 3, ntests);
  SET_VECTOR_ELT(learning, 4, algo);
  SET_VECTOR_ELT(learning, 5, args);

  /* allocate and initialize the main structure. */
  PROTECT(res = allocVector(VECSXP, 3));
  setAttrib(res, R_NamesSymbol, names2);
  setAttrib(res, R_ClassSymbol, class);
  SET_VECTOR_ELT(res, 0, learning);
  SET_VECTOR_ELT(res, 1, cached);
  SET_VECTOR_ELT(res, 2, arcs);

  UNPROTECT(8);
  return res;

}/*BN_BASE_STRUCTURE*/

static int ic_logic(int *amat, SEXP nodes, int *nnodes, int *arc, int *work,
    int *degree, double *max, int *in_degree, double *max_in, int *out_degree,
    double *max_out, int *cozman, int *debuglevel) {

int path = 0;

  /* sample an arc (that is, two nodes different from each other). */
  SampleNoReplace(2, *nnodes, arc, work);

  if (amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] == 1) {

    /* if the arc (i, j) exists in the actual graph, delete the arc,
     * provided that the underlying graph remains connected. */
    if (*debuglevel > 0)
      Rprintf("  > arc %s -> %s is present.\n",
        NODE(arc[0] - 1), NODE(arc[1] - 1));

    /* skip the check on connectedness for Melancon's algorithm, it's not
     * required.  */
    if (*cozman == 1)  {

      /* if there is a(n undirected) path in the (underlying undirected) graph
       * from arc[0] to arc[1] other than arc[0] -> arc[1], the graph is still
       * connected. */
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 0;
      path = c_has_path(arc[0] - 1, arc[1] - 1, amat, *nnodes, nodes,
               TRUE, FALSE, FALSE);
      amat[CMC(arc[0] - 1, arc[1] - 1, *nnodes)] = 1;

    }/*THEN*/
    else {

      path = 1;

    }/*ELSE*/

    if (path) {

      if (*debuglevel > 0)
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

      if (*debuglevel > 0)
        Rprintf("  @ not removing arc %s -> %s (graph not connected).\n",
          NODE(arc[0] - 1), NODE(arc[1] - 1));

      return 0;

    }/*ELSE*/

  }/*THEN*/
  else {

    /* add the arc, provided that the underlying graph remains acyclic. */
    if (*debuglevel > 0)
      Rprintf("  > arc %s -> %s is not present.\n", NODE(arc[0] - 1), NODE(arc[1] - 1));

    /* do not add the arc if this violates the constraints on the degrees of the
     * nodes it's incident on. */
    if ((degree[arc[0] - 1] >= *max) || (degree[arc[1] - 1] >= *max) ||
        (out_degree[arc[0] - 1] >= *max_out) || (in_degree[arc[1]] >= *max_in)) {

      if (*debuglevel > 0) {

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
    path = c_has_path(arc[1] - 1, arc[0] - 1, amat, *nnodes,
                     nodes, FALSE, FALSE, FALSE);

    if (!path) {

      if (*debuglevel > 0)
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

      if (*debuglevel > 0)
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
  SETCAR(t, install("modelstring"));
  t = CDR(t);
  /* second slot, the bayesian network (the only argument). */
  SETCAR(t,  bn);
  /* evaluate ... */
  PROTECT(t = eval(s, R_GlobalEnv));
  /* ... and print the result. */
  Rprintf("  > model string is:\n%s\n", CHAR(STRING_ELT(t, 0)));
  UNPROTECT(2);

}/*PRINT_MODELSTRING*/

