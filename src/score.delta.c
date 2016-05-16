#include "include/rcore.h"
#include "include/globals.h"
#include "include/scores.h"

static SEXP score_delta_helper(SEXP net, SEXP arc, SEXP operator, int children,
    int both) {

int i = 0, k = 0;
char *from = (char *) CHAR(STRING_ELT(arc, 0));
char *to = (char *) CHAR(STRING_ELT(arc, 1));
char *op = (char *) CHAR(STRING_ELT(operator, 0));
SEXP fake, nodes, parents_from, parents_to, children_from;
SEXP start, cur_n, cur_p, cur_c, temp, name_nodes, name_cached;

  /* allocate the return value. */
  PROTECT(fake = allocVector(VECSXP, 1));

  /* allocate and initialize the names of the elements of the fake network. */
  PROTECT(name_nodes = mkString("nodes"));
  PROTECT(name_cached = allocVector(STRSXP, 1 + children));
  SET_STRING_ELT(name_cached, 0, mkChar("parents"));
  if (children)
    SET_STRING_ELT(name_cached, 1, mkChar("children"));

  start = getListElement(net, "nodes");

  if (!strcmp(op, "set")) {

    /* adding the arc to the graph. */
    PROTECT(nodes = allocVector(VECSXP, 1 + both));
    PROTECT(temp = allocVector(VECSXP, 1 + children));

    cur_n = getListElement(start, to);
    cur_p = getListElement(cur_n, "parents");

    PROTECT(parents_to = allocVector(STRSXP, length(cur_p) + 1));

    /* the new parents are the old ones plus the other node incident on
     * the arc. */
    for (i = 0; i < length(cur_p); i++)
      SET_STRING_ELT(parents_to, i, STRING_ELT(cur_p, i));
    SET_STRING_ELT(parents_to, length(cur_p), STRING_ELT(arc, 0));

    SET_VECTOR_ELT(temp, 0, parents_to);
    if (children)
      SET_VECTOR_ELT(temp, 1, getListElement(cur_n, "children"));

    if (both) {

      /* set the names in the structure. */
      setAttrib(nodes, R_NamesSymbol, arc);
      setAttrib(temp, R_NamesSymbol, name_cached);
      /* assign the second node. */
      SET_VECTOR_ELT(nodes, 1, duplicate(temp));

      /* work on the second node. */
      cur_n = getListElement(start, from);
      cur_c = getListElement(cur_n, "children");

      PROTECT(children_from = allocVector(STRSXP, length(cur_c) + 1));

      /* the new children are the old ones plus the other node incident on
       * the arc. */
      for (i = 0; i < length(cur_c); i++)
        SET_STRING_ELT(children_from, i, STRING_ELT(cur_c, i));
      SET_STRING_ELT(children_from, length(cur_c), STRING_ELT(arc, 1));

      SET_VECTOR_ELT(temp, 0, getListElement(cur_n, "parents"));
      SET_VECTOR_ELT(temp, 1, children_from);
      /* assign the first node. */
      SET_VECTOR_ELT(nodes, 0, temp);

      UNPROTECT(1);

    }/*THEN*/
    else {

      /* set the names in the structure. */
      setAttrib(nodes, R_NamesSymbol, mkString(to));
      setAttrib(temp, R_NamesSymbol, name_cached);
      /* assign the only node. */
      SET_VECTOR_ELT(nodes, 0, temp);

    }/*ELSE*/

    UNPROTECT(3);

  }/*THEN*/
  else if (!strcmp(op, "drop")) {

    /* dropping the arc from the graph. */
    PROTECT(nodes = allocVector(VECSXP, 1 + both));
    PROTECT(temp = allocVector(VECSXP, 1 + children));

    cur_n = getListElement(start, to);
    cur_p = getListElement(cur_n, "parents");

    PROTECT(parents_to = allocVector(STRSXP, length(cur_p) - 1));

    /* the new parents are the old ones except the other node incident on
     * the arc. */
    for (i = 0, k = 0; i < length(cur_p); i++)
      if (strcmp(CHAR(STRING_ELT(cur_p, i)), from) != 0)
        SET_STRING_ELT(parents_to, k++, STRING_ELT(cur_p, i));

    SET_VECTOR_ELT(temp, 0, parents_to);
    if (children)
      SET_VECTOR_ELT(temp, 1, getListElement(cur_n, "children"));

    if (both) {

      /* set the names in the structure. */
      setAttrib(nodes, R_NamesSymbol, arc);
      setAttrib(temp, R_NamesSymbol, name_cached);
      /* assign the second node. */
      SET_VECTOR_ELT(nodes, 1, duplicate(temp));

      /* work on the second node. */
      cur_n = getListElement(start, from);
      cur_c = getListElement(cur_n, "children");

      PROTECT(children_from = allocVector(STRSXP, length(cur_c) - 1));

      /* the new children are the old ones except the other node incident on
       * the arc. */
      for (i = 0, k = 0; i < length(cur_c); i++)
        if (strcmp(CHAR(STRING_ELT(cur_c, i)), to) != 0)
          SET_STRING_ELT(children_from, k++, STRING_ELT(cur_c, i));

      SET_VECTOR_ELT(temp, 0, getListElement(cur_n, "parents"));
      SET_VECTOR_ELT(temp, 1, children_from);
      /* assign the first node. */
      SET_VECTOR_ELT(nodes, 0, temp);

      UNPROTECT(1);

    }/*THEN*/
    else {

      /* set the names in the structure. */
      setAttrib(nodes, R_NamesSymbol, mkString(to));
      setAttrib(temp, R_NamesSymbol, name_cached);
      /* assign the only node. */
      SET_VECTOR_ELT(nodes, 0, temp);

    }/*ELSE*/

    UNPROTECT(3);

  }/*THEN*/
  else {

    /* reversing the arc in the graph. */
    PROTECT(nodes = allocVector(VECSXP, 2));
    PROTECT(temp = allocVector(VECSXP, 1 + children));
    setAttrib(nodes, R_NamesSymbol, arc);
    setAttrib(temp, R_NamesSymbol, name_cached);

    /* add "to" to the parents of "from". */
    cur_n = getListElement(start, from);
    cur_p = getListElement(cur_n, "parents");

    PROTECT(parents_from = allocVector(STRSXP, length(cur_p) + 1));

    for (i = 0; i < length(cur_p); i++)
      SET_STRING_ELT(parents_from, i, STRING_ELT(cur_p, i));
    SET_STRING_ELT(parents_from, length(cur_p), STRING_ELT(arc, 1));

    SET_VECTOR_ELT(temp, 0, parents_from);
    if (children)
      SET_VECTOR_ELT(temp, 1, getListElement(cur_n, "children"));
    SET_VECTOR_ELT(nodes, 0, duplicate(temp));

    /* remove "from" from the parents of "to". */
    cur_n = getListElement(start, to);
    cur_p = getListElement(cur_n, "parents");

    PROTECT(parents_to = allocVector(STRSXP, length(cur_p) - 1));

    for (i = 0, k = 0; i < length(cur_p); i++)
      if (strcmp(CHAR(STRING_ELT(cur_p, i)), from) != 0)
        SET_STRING_ELT(parents_to, k++, STRING_ELT(cur_p, i));

    SET_VECTOR_ELT(temp, 0, parents_to);
    if (children)
      SET_VECTOR_ELT(temp, 1, getListElement(cur_n, "children"));
    SET_VECTOR_ELT(nodes, 1, temp);

    UNPROTECT(4);

  }/*ELSE*/

  /* save the fabricated nodes' structures in the return value. */
  SET_VECTOR_ELT(fake, 0, nodes);
  setAttrib(fake, R_NamesSymbol, name_nodes);

  UNPROTECT(3);

  return fake;

}/*SCORE_DELTA_HELPER*/

static double robust_score_difference(double old1, double old2, double new1,
    double new2) {

double delta = 0, new_sum = 0, old_sum = 0;

  /* sum old and new score components. */
  new_sum = new1 + new2;
  old_sum = old1 + old2;

  /* compare the network scores, minus numeric tolerance for better score
   * equivalence detection. */
  if (fabs(new_sum - old_sum) < MACHINE_TOL)
    delta = 0;
  else
    delta = new_sum - old_sum;

  /* catch: the difference between two -Inf scores must be -Inf and not NaN,
   * so that we can safely test it and reject the change; and if the old score
   * is -Inf just look at the new score. */
  if (old_sum == R_NegInf)
    delta = (new_sum != R_NegInf) ? new_sum : R_NegInf;

  return delta;

}/*ROBUST_DIFFERENCE*/

SEXP score_delta_decomposable(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra, int chld) {

int *t = NULL;
double diff = 0, *new = NULL, *old = NULL;
SEXP delta, fake, new_score, try, to_update;

  /* create the fake network with the updated structure. */
  PROTECT(fake = score_delta_helper(network, arc, op, chld, FALSE));
  /* find out which nodes to update from the fake structure. */
  to_update = getAttrib(getListElement(fake, "nodes"), R_NamesSymbol);

  /* compute the updated scores with the fake newtork. */
  PROTECT(new_score = allocVector(REALSXP, length(to_update)));
  new = REAL(new_score);
  c_per_node_score(fake, data, score, to_update, extra, FALSE, new);
  /* update the test counter. */
  test_counter += length(new_score);

  /* get the corresponding components from the old score. */
  PROTECT(try = match(getAttrib(reference_score, R_NamesSymbol), to_update, 0));
  t = INTEGER(try);
  old = REAL(reference_score);

  if (length(new_score) == 1)
    diff = robust_score_difference(old[t[0] - 1], 0, new[0], 0);
  else
    diff = robust_score_difference(old[t[0] - 1], old[t[1] - 1], new[0], new[1]);

  UNPROTECT(1);

  /* build the return value. */
  PROTECT(delta = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(delta, 0, ScalarLogical(diff > 0));
  SET_VECTOR_ELT(delta, 1, ScalarReal(diff));
  SET_VECTOR_ELT(delta, 2, new_score);
  setAttrib(delta, R_NamesSymbol, mkStringVec(3, "bool", "delta", "updates"));

  UNPROTECT(3);

  return delta;

}/*SCORE_DELTA_DECOMPOSABLE*/

SEXP score_delta_cs(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra) {

int *t = NULL;
const char *o = CHAR(STRING_ELT(op, 0));
double retval = 0, *new = 0, *old = NULL, new_prior = 0, old_prior = 0;
SEXP fake, prior, cs_prior, order, try, delta, new_score;
SEXP fake_nodes, real_nodes, target;

  /* get the prior specification. */
  prior = getListElement(extra, "prior");
  cs_prior = getListElement(extra, "beta");
  order = getAttrib(cs_prior, BN_NodesSymbol);

  PROTECT(try = match(order, arc, 0));
  t = INTEGER(try);

  if ((t[0] > t[1]) || (strcmp(o, "reverse") == 0)) {

    UNPROTECT(1);

    return score_delta_decomposable(arc, network, data, score, score_delta,
             reference_score, op, extra, TRUE);

  }/*THEN*/
  else {

    /* create a fake network with the updated structure, including children. */
    PROTECT(fake = score_delta_helper(network, arc, op, TRUE, TRUE));
    /* find out which nodes to update from the fake structure. */
    fake_nodes = getListElement(fake, "nodes");
    PROTECT(target = allocVector(STRSXP, 1));
    SET_STRING_ELT(target, 0, STRING_ELT(arc, 1));
    /* compute the one score component that needs recomputing. */
    PROTECT(new_score = allocVector(REALSXP, 2));
    new = REAL(new_score);
    c_per_node_score(fake, data, score, target, extra, FALSE, new + 1);
    /* update the test counter. */
    test_counter++;

    /* update the prior of the other one. */
    SET_STRING_ELT(target, 0, STRING_ELT(arc, 0));

    fake_nodes = VECTOR_ELT(fake_nodes, 0);
    new_prior = graph_prior_prob(prior, target, fake_nodes, cs_prior, FALSE);
    real_nodes = VECTOR_ELT(getListElement(network, "nodes"), t[0] - 1);
    old_prior = graph_prior_prob(prior, target, real_nodes, cs_prior, FALSE);

    old = REAL(reference_score);
    new[0] = old[t[0] - 1] - old_prior + new_prior;

    /* compute the difference. */
    if (length(new_score) == 1)
      retval = robust_score_difference(old[t[0] - 1], 0, new[0], 0);
    else
      retval = robust_score_difference(old[t[0] - 1], old[t[1] - 1], new[0], new[1]);

    /* build the return value. */
    PROTECT(delta = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(delta, 0, ScalarLogical(retval > 0));
    SET_VECTOR_ELT(delta, 1, ScalarReal(retval));
    SET_VECTOR_ELT(delta, 2, new_score);
    setAttrib(delta, R_NamesSymbol, mkStringVec(3, "bool", "delta", "updates"));

    UNPROTECT(5);

    return delta;

  }/*ELSE*/

}/*SCORE_DELTA_CS*/

SEXP score_delta_monolithic(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra) {

const char *s = CHAR(STRING_ELT(score, 0));

  if ((strcmp(s, "bde") == 0) || (strcmp(s, "bge") == 0) ||
      (strcmp(s, "bds") == 0) || (strcmp(s, "mbde") == 0)) {

    const char *prior = CHAR(STRING_ELT(getListElement(extra, "prior"), 0));

    if (strcmp(prior, "cs") == 0)
      return score_delta_cs(arc, network, data, score, score_delta,
               reference_score, op, extra);

  }/*THEN*/

  return R_NilValue;

}/*SCORE_DELTA_MONOLITHIC*/

SEXP score_delta(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra, SEXP decomposable) {

  if (isTRUE(decomposable)) {

    return score_delta_decomposable(arc, network, data, score, score_delta,
             reference_score, op, extra, FALSE);

  }/*THEN*/
  else {

    return score_delta_monolithic(arc, network, data, score, score_delta,
             reference_score, op, extra);

  }/*ELSE*/

}/*SCORE_DELTA*/
