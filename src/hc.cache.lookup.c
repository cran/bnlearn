#include "common.h"

SEXP hc_cache_fill(SEXP nodes, SEXP data, SEXP network, SEXP score, SEXP extra,
    SEXP reference, SEXP equivalence, SEXP updated, SEXP env, SEXP amat,
    SEXP cache, SEXP blmat, SEXP debug) {

int *colsum = NULL, nnodes = LENGTH(nodes), lupd = LENGTH(updated);
int *a = NULL, *upd = NULL, *b = NULL, *debuglevel = NULL;
int i = 0, j = 0, k = 0;
double *cache_value = NULL;
SEXP arc, callenv, params, delta, op, dummy, temp;

  /* dereference the debug parameter. */
  debuglevel = LOGICAL(debug);

  /* save a pointer to the adjacency matrix, the blacklist and the
   * updated nodes. */
  a = INTEGER(amat);
  b = INTEGER(blmat);
  upd = INTEGER(updated);

  /* if there are no nodes to update, return. */
  if (lupd == 0) return cache;

  /* set up row and column total to check for score equivalence;
   * zero means no parent nodes. */
  if (isTRUE(equivalence)) {

    colsum = alloc1dcont(nnodes);

    for (i = 0; i < nnodes; i++)
      for (j = 0; j < nnodes; j++)
        colsum[j] += a[CMC(i, j, nnodes)];

  }/*THEN*/

  /* allocate and initialize the cache. */
  cache_value = REAL(cache);

  /* allocate a two-slot character vector. */
  PROTECT(arc = allocVector(STRSXP, 2));

  /* allocate the call environment fro the call to score.delta. */
  PROTECT(params = callenv = allocList(10));
  SET_TYPEOF(callenv, LANGSXP);

  /* allocate and initialize the fake score delta. */
  PROTECT(delta = allocVector(REALSXP, 1));
  NUM(delta) = 0;

  /* allocate and initialize the score.delta() operator. */
  PROTECT(op = allocVector(STRSXP, 1));
  SET_STRING_ELT(op, 0, mkChar("set"));

  /* allocate and initialize the debug parameter. */
  PROTECT(dummy = allocVector(LGLSXP, 1));
  LOGICAL(dummy)[0] = FALSE;

  /* set up the call to score.delta(). */
  SETCAR(params, install("score.delta"));
  params = CDR(params);
  SETCAR(params, arc);
  params = CDR(params);
  SETCAR(params, network);
  params = CDR(params);
  SETCAR(params, data);
  params = CDR(params);
  SETCAR(params, score);
  params = CDR(params);
  SETCAR(params, delta);
  params = CDR(params);
  SETCAR(params, reference);
  params = CDR(params);
  SETCAR(params, op);
  params = CDR(params);
  SETCAR(params, extra);
  params = CDR(params);
  SETCAR(params, dummy);

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

       /* incident nodes must be different from each other. */
       if (i == j) continue;

       /* if only one or two nodes' caches need updating, skip the rest. */
       for (k = 0; k < lupd; k++)
         if (upd[k] == j)
           goto there;

       continue;

there:

       /* no need to compute the score delta for blacklisted arcs. */
       if (b[CMC(i, j, nnodes)] == 1)
         continue;

       /* use score equivalence if possible to check only one orientation. */
       if (isTRUE(equivalence)) {

         /* if the following conditions are met, look up the score delta of
          * the reverse of the current arc:
          *   1) that score delta has already been computed.
          *   2) both incident nodes have no parent, so the arc is really
          *      score equivalent (no v-structures).
          *   3) the reversed arc has not been blacklisted, as the score delta
          *      is not computed in this case. */
         if ((i > j) && (colsum[i] + colsum[j] == 0) && (b[CMC(j, i, nnodes)] == 0)) {

           cache_value[CMC(i, j, nnodes)] = cache_value[CMC(j, i, nnodes)];
           continue;

         }/*THEN*/

       }/*THEN*/

       /* save the nodes incident on the arc. */
       SET_STRING_ELT(arc, 0, STRING_ELT(nodes, i));
       SET_STRING_ELT(arc, 1, STRING_ELT(nodes, j));

       /* add the arc to the call to score.delta().  */
       params = CDR(callenv);
       SETCAR(params, arc);

       /* if the arc is not present in the graph it should be added;
        * otherwise it should be removed. */
       for (k = 0; k < 6; k++) params = CDR(params);

       if (a[CMC(i, j, nnodes)] == 0)
         SET_STRING_ELT(op, 0, mkChar("set"));
       else
         SET_STRING_ELT(op, 0, mkChar("drop"));

       SETCAR(params, op);

       /* evaluate the call to score.delta() for the arc. */
       PROTECT(temp = eval(callenv, env));
       cache_value[CMC(i, j, nnodes)] = NUM(VECTOR_ELT(temp, 1));
       UNPROTECT(1);

       if (*debuglevel > 0)
         Rprintf("* caching score delta for arc %s -> %s (%lf).\n",
           CHAR(STRING_ELT(nodes, i)), CHAR(STRING_ELT(nodes, j)),
            cache_value[CMC(i, j, nnodes)]);

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(5);
  return cache;

}/*HC_CACHE_FILL*/

/* a single step of the optimized hill climbing (one arc addition/removal/reversal). */
SEXP hc_opt_step(SEXP amat, SEXP nodes, SEXP added, SEXP cache, SEXP reference,
    SEXP wlmat, SEXP blmat, SEXP debug) {

int nnodes = LENGTH(nodes), i = 0, j = 0;
int *am = NULL, *ad = NULL, *w = NULL, *b = NULL, *debuglevel = NULL;
int counter = 0, update = 1, from = 0, to = 0;
double *cache_value = NULL, temp = 0, max = 0, tol = MACHINE_TOL;
SEXP false, bestop, names;

  /* allocate and initialize the return value (use FALSE as a canary value). */
  PROTECT(bestop = allocVector(VECSXP, 3));
  PROTECT(names = allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("op"));
  SET_STRING_ELT(names, 1, mkChar("from"));
  SET_STRING_ELT(names, 2, mkChar("to"));
  setAttrib(bestop, R_NamesSymbol, names);

  /* allocate and initialize a dummy FALSE object. */
  PROTECT(false = allocVector(LGLSXP, 1));
  LOGICAL(false)[0] = FALSE;
  SET_VECTOR_ELT(bestop, 0, false);

  /* save pointers to the numeric/integer matrices. */
  cache_value = REAL(cache);
  ad = INTEGER(added);
  am = INTEGER(amat);
  w = INTEGER(wlmat);
  b = INTEGER(blmat);

  /* dereference the debug parameter. */
  debuglevel = LOGICAL(debug);

  if (*debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0; i < nnodes * nnodes; i++)
       counter += ad[i];

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to add one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (ad[CMC(i, j, nnodes)] == 0)
        continue;

      /* retrieve the score delta from the cache. */
      temp = cache_value[CMC(i, j, nnodes)];

      if (*debuglevel > 0) {

        Rprintf("  > trying to add %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp);

      }/*THEN*/

      /* this score delta is the best one at the moment, so add the arc if it
       * does not introduce cycles in the graph. */
      if (temp - max > tol) {

        if (c_has_path(j, i, am, nnodes, nodes, FALSE, FALSE, FALSE)) {

          if (*debuglevel > 0)
            Rprintf("    > not adding, introduces cycles in the graph.\n");

          continue;

        }/*THEN*/

        if (*debuglevel > 0)
          Rprintf("    @ adding %s -> %s.\n", NODE(i), NODE(j));

        /* update the return value. */
        bestop_update(bestop, "set", NODE(i), NODE(j));
        /* store the node indices to update the reference scores. */
        from = i;
        to = j;

        /* update the threshold score delta. */
        max = temp;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  if (*debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0, counter = 0; i < nnodes * nnodes; i++)
       counter += am[i] * (1 - w[i]);

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to remove one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (am[CMC(i, j, nnodes)] == 0)
        continue;

      /* whitelisted arcs are not to be removed, ever. */
      if (w[CMC(i, j, nnodes)] == 1)
        continue;

      /* retrieve the score delta from the cache. */
      temp = cache_value[CMC(i, j, nnodes)];

      if (*debuglevel > 0) {

        Rprintf("  > trying to remove %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp);

      }/*THEN*/

      if (temp - max > tol) {

        if (*debuglevel > 0)
          Rprintf("    @ removing %s -> %s.\n", NODE(i), NODE(j));

        /* update the return value. */
        bestop_update(bestop, "drop", NODE(i), NODE(j));
        /* store the node indices to update the reference scores. */
        from = i;
        to = j;

        /* update the threshold score delta. */
        max = temp;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  if (*debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0, counter = 0; i < nnodes; i++)
       for (j = 0; j < nnodes; j++)
         counter += am[CMC(i, j, nnodes)] * (1 - b[CMC(j, i, nnodes)]);

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to reverse one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (am[CMC(i, j, nnodes)] == 0)
        continue;

      /* don't reverse an arc if the one in the opposite direction is
       * blacklisted, ever. */
      if (b[CMC(j, i, nnodes)] == 1)
        continue;

      /* retrieve the score delta from the cache. */
      temp = cache_value[CMC(i, j, nnodes)] + cache_value[CMC(j, i, nnodes)];

      if (*debuglevel > 0) {

        Rprintf("  > trying to reverse %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp);

      }/*THEN*/

      if (temp - max > tol) {

        if (c_has_path(i, j, am, nnodes, nodes, FALSE, TRUE, FALSE)) {

          if (*debuglevel > 0)
            Rprintf("    > not reversing, introduces cycles in the graph.\n");

          continue;

        }/*THEN*/

        if (*debuglevel > 0)
          Rprintf("    @ reversing %s -> %s.\n", NODE(i), NODE(j));

        /* update the return value. */
        bestop_update(bestop, "reverse", NODE(i), NODE(j));
        /* store the node indices to update the reference scores. */
        from = i;
        to = j;
        /* both nodes' reference scores must be updated. */
        update = 2;

        /* update the threshold score delta. */
        max = temp;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  /* update the reference scores. */
  REAL(reference)[to] += cache_value[CMC(from, to, nnodes)];
  if (update == 2)
    REAL(reference)[from] += cache_value[CMC(to, from, nnodes)];

  UNPROTECT(3);

  return bestop;

}/*HC_OPT_STEP*/

void bestop_update(SEXP bestop, char *op, const char *from, const char *to) {

SEXP temp;

  PROTECT(temp = allocVector(STRSXP, 1));

  SET_STRING_ELT(temp, 0, mkChar(op));
  SET_VECTOR_ELT(bestop, 0, duplicate(temp));

  SET_STRING_ELT(temp, 0, mkChar(from));
  SET_VECTOR_ELT(bestop, 1, duplicate(temp));

  SET_STRING_ELT(temp, 0, mkChar(to));
  SET_VECTOR_ELT(bestop, 2, duplicate(temp));

  UNPROTECT(1);

}/*BESTOP_UPDATE*/

SEXP hc_to_be_added(SEXP arcs, SEXP blacklist, SEXP whitelist, SEXP nodes,
    SEXP convert) {

int i = 0, j = 0, narcs = 0, dims = LENGTH(nodes);
int *a = NULL, *coords = NULL;
short int duplicated = 0;
SEXP try, result = R_NilValue, result2;

  /* transform the arc set into an adjacency matrix, if it's not one already. */
  if (isInteger(arcs)) {

    if ((duplicated = NAMED(arcs)) > 0)
      PROTECT(result = duplicate(arcs));

  }/*THEN*/
  else {

    PROTECT(result = arcs2amat(arcs, nodes));

  }/*ELSE*/

  /* dereference the adjacency matrix once and for all. */
  a = INTEGER(result);

  /* flip all the nondiagonal cells. */
  for (j = 0; j < dims; j++) {

    for (i = 0; i < dims; i++) {

      /* diagonal elements are always equal to zero, skip them. */
      if (i == j)
        continue;

      a[CMC(i, j, dims)] = 1 - a[CMC(i, j, dims)];

    }/*FOR*/

  }/*FOR*/

  /* if an arc cannot be added in one direction, it cannot be added in the
   * opposite one either; flip them off the adjacency matrix. */
  for (j = 0; j < dims; j++)
    for (i = j + 1; i < dims; i++)
      a[CMC(j, i, dims)] = a[CMC(i, j, dims)] = a[CMC(i, j, dims)] * a[CMC(j, i, dims)];

#define FLIP_FROM_LIST(list, value) \
  if (!isNull(list)) { \
    if (!isInteger(list)) { \
      PROTECT(try = match(nodes, list, 0)); \
      coords = INTEGER(try); \
      narcs = LENGTH(try)/2; \
      for (i = 0; i < narcs; i++)  \
        a[CMC(coords[i] - 1, coords[i + narcs] - 1, dims)] = value; \
      UNPROTECT(1); \
    } \
    else { \
      coords = INTEGER(list); \
      for (i = 0; i < dims * dims; i ++) \
        if (coords[i] == 1) \
          a[i] = value; \
    } \
  }

  /* now the blacklist gets involved. */
  FLIP_FROM_LIST(blacklist, 0);
  /* and, last but not least, the whitelist gets involved. */
  FLIP_FROM_LIST(whitelist, 1);

  /* return either the adjacency matrix or the arc set. */
  if (isTRUE(convert)) {

    PROTECT(result2 = amat2arcs(result, nodes));

    if ((duplicated > 0) || !isInteger(arcs))
      UNPROTECT(2);
    else
      UNPROTECT(1);
    return result2;

  }/*THEN*/
  else {

    if ((duplicated > 0) || !isInteger(arcs))
      UNPROTECT(1);
    return result;

  }/*ELSE*/

}/*HC_TO_BE_ADDED*/

