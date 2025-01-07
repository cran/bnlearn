#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/uppertriangular.h"
#include "../core/sort.h"
#include "../minimal/data.frame.h"
#include "../minimal/strings.h"
#include "../minimal/unique.h"
#include "../minimal/common.h"
#include "../include/globals.h"
#include "../include/graph.h"
#include "scores.h"

double castelo_prior(SEXP beta, SEXP target, SEXP parents, SEXP children,
    bool debugging);
double marginal_prior(SEXP beta, SEXP target, SEXP parents, SEXP children,
    SEXP node_cache, SEXP nodes, bool debugging);

double graph_prior_prob(SEXP prior, SEXP target, SEXP beta, SEXP cache,
    bool debugging) {

double *b = NULL, prob = 0;
SEXP nodes, parents, children, target_cache;

  /* check which prior should be computed, and use the uniform one
   * if none is specified.*/
  if (prior ==  R_NilValue)
    return 0;

  /* extract the cached information about the target node. */
  target_cache = getListElement(cache, (char *)CHAR(STRING_ELT(target, 0)));

  switch(gprior_to_enum(CHAR(STRING_ELT(prior, 0)))) {

    /* constant prior, log(1) = 0 for backward compatibility. */
    case UNIFORM:
      prob = 0;
      break;

    /* variable selection prior, each arc has independent beta probability
     * for inclusion. */
    case VSP:
      parents = getListElement(target_cache, "parents");
      b = REAL(beta);
      prob = length(parents) * log(*b / (1 - *b));
      break;

    /* completed prior from Castelo and Siebes. */
    case CS:
      parents = getListElement(target_cache, "parents");
      children = getListElement(target_cache, "children");

      if (beta == R_NilValue)
        prob = 0;
      else
        prob = castelo_prior(beta, target, parents, children, debugging);
      break;

    /* marignal uniform prior. */
    case MU:
      parents = getListElement(target_cache, "parents");
      children = getListElement(target_cache, "children");
      nodes = getAttrib(beta, BN_NodesSymbol);
      prob = marginal_prior(beta, target, parents, children, cache, nodes,
               debugging);
      break;

    default:
      error("unknown graph prior.");

  }/*SWITCH*/

  return prob;

}/*GRAPH_PRIOR_PROB*/

#define PARENT 1
#define CHILD  2

double castelo_prior(SEXP beta, SEXP target, SEXP parents, SEXP children,
    bool debugging) {

int i = 0, k = 0, t = 0, nnodes = 0, cur_arc = 0;
int nbeta = length(VECTOR_ELT(beta, 0));
int *temp = NULL, *aid = INTEGER(VECTOR_ELT(beta, 2));
double prior = 0, result = 0;
double *bkwd = REAL(VECTOR_ELT(beta, 4)), *fwd = REAL(VECTOR_ELT(beta, 3));
short int *adjacent = NULL;
SEXP nodes, try;

  /* get the node labels. */
  nodes = getAttrib(beta, BN_NodesSymbol);
  nnodes = length(nodes);

  /* match the target node. */
  PROTECT(try = match(nodes, target, 0));
  t = INT(try);
  UNPROTECT(1);

  /* find out which nodes are parents and which nodes are children. */
  adjacent = Calloc1D(nnodes, sizeof(short int));

  PROTECT(try = match(nodes, parents, 0));
  temp = INTEGER(try);
  for (i = 0; i < length(try); i++)
    adjacent[temp[i] - 1] = PARENT;
  UNPROTECT(1);

  PROTECT(try = match(nodes, children, 0));
  temp = INTEGER(try);
  for (i = 0; i < length(try); i++)
    adjacent[temp[i] - 1] = CHILD;
  UNPROTECT(1);

  /* prior probabilities table lookup. */
  for (i = t + 1; i <= nnodes; i++) {

    /* compute the arc id. */
    cur_arc = UPTRI3(t, i, nnodes);

    /* look up the prior probability. */
    for (/*k,*/ prior = ((double)1/3); k < nbeta; k++) {

      /* arcs are ordered, so we can stop early in the lookup. */
      if (aid[k] > cur_arc)
        break;

      if (aid[k] == cur_arc) {

        switch(adjacent[i - 1]) {

          case PARENT:
            prior = bkwd[k];
            break;
          case CHILD:
            prior = fwd[k];
            break;
          default:
            prior = fmax2(0, 1 - bkwd[k] - fwd[k]);

        }/*SWITCH*/

        break;

      }/*THEN*/

    }/*FOR*/

    if (debugging) {

      switch(adjacent[i - 1]) {

        case PARENT:
          Rprintf("  > found arc %s -> %s, prior pobability is %lf.\n",
            NODE(i - 1), NODE(t - 1), prior);
          break;
        case CHILD:
          Rprintf("  > found arc %s -> %s, prior probability is %lf.\n",
            NODE(t - 1), NODE(i - 1), prior);
          break;
        default:
          Rprintf("  > no arc between %s and %s, prior probability is %lf.\n",
            NODE(t - 1), NODE(i - 1), prior);

      }/*SWITCH*/

    }/*THEN*/

    /* move to log-scale and divide by the non-informative log(1/3), so that
     * the contribution of each arc whose prior has not been specified by the
     * user is zero; overflow is likely otherwise. */
    result += log(prior / ((double)1/3));

  }/*FOR*/

  Free1D(adjacent);

  return result;

}/*CASTELO_PRIOR*/

/* complete a prior as per Castelo & Siebes. */
SEXP castelo_completion(SEXP prior, SEXP nodes, SEXP learning)  {

int i = 0, k = 0, cur = 0, narcs1 = 0, narcs2 = 0, nnodes = length(nodes);
int *m1 = NULL, *m2 = NULL, *und = NULL, *aid = NULL, *poset = NULL, *id = NULL;
double *d1 = NULL, *d2 = NULL, *p = NULL, tol = MACHINE_TOL;
SEXP df, arc_id, undirected, a1, a2, match1, match2, prob;
SEXP result, from, to, nid, dir1, dir2;

  /* compute numeric IDs for the arcs. */
  a1 = VECTOR_ELT(prior, 0);
  a2 = VECTOR_ELT(prior, 1);
  narcs1 = length(a1);
  PROTECT(match1 = match(nodes, a1, 0));
  PROTECT(match2 = match(nodes, a2, 0));
  m1 = INTEGER(match1);
  m2 = INTEGER(match2);
  PROTECT(arc_id = allocVector(INTSXP, narcs1));
  aid = INTEGER(arc_id);

  c_arc_hash(narcs1, nnodes, m1, m2, aid, NULL, FALSE);

  /* duplicates correspond to undirected arcs. */
  PROTECT(undirected = dupe(arc_id));
  und = INTEGER(undirected);

  /* extract the components from the prior. */
  prob = VECTOR_ELT(prior, 2);
  p = REAL(prob);

  /* count output arcs. */
  for (i = 0; i < narcs1; i++)
    narcs2 += 2 - und[i];
  narcs2 /= 2;

  /* allocate the columns of the return value. */
  PROTECT(from = allocVector(STRSXP, narcs2));
  PROTECT(to = allocVector(STRSXP, narcs2));
  PROTECT(nid = allocVector(INTSXP, narcs2));
  id = INTEGER(nid);
  PROTECT(dir1 = allocVector(REALSXP, narcs2));
  d1 = REAL(dir1);
  PROTECT(dir2 = allocVector(REALSXP, narcs2));
  d2 = REAL(dir2);

  /* sort the strength coefficients. */
  poset = Calloc1D(narcs1, sizeof(int));
  for (k = 0; k < narcs1; k++)
    poset[k] = k;
  i_sort(aid, poset, narcs1);

  for (i = 0, k = 0; i < narcs1; i++) {

    cur = poset[i];

#define ASSIGN(A1, A2, D1, D2) \
  SET_STRING_ELT(from,  k, STRING_ELT(A1, cur)); \
  SET_STRING_ELT(to,  k, STRING_ELT(A2, cur)); \
  id[k] = aid[i]; \
  D1[k] = p[cur]; \
  if ((und[cur] == TRUE) && (i < narcs1 - 1)) \
    D2[k] = p[poset[++i]]; \
  else \
    D2[k] = (1 - D1[k])/2;

    /* copy the node labels. */
    if (m1[cur] < m2[cur]) {

      ASSIGN(a1, a2, d1, d2);

    }/*THEN*/
    else {

      ASSIGN(a2, a1, d2, d1);

    }/*ELSE*/

    /* check the probabilities do not exceed 1; fail only for large errors. */
    if (d1[k] + d2[k] > 1) {

      if (d1[k] + d2[k] < 1.01) {

        d1[k] = d1[k] / (d1[k] + d2[k]);
        d2[k] = d2[k] / (d1[k] + d2[k]);

      }/*THEN*/
      else {

        UNPROTECT(9);

        error("the probabilities for arc %s -> %s sum to %lf.",
          CHAR(STRING_ELT(from, k)), CHAR(STRING_ELT(to, k)), d1[k] + d2[k]);

      }/*ELSE*/

    }/*THEN*/

    /* shrink the probabilities away from 0 and 1, structure learning otherwise
     * fails when starting from the empty graph and gets stuck very easily when
     * starting from a non-empty graph (and in general). */
    if (isTRUE(learning)) {

      if ((d1[k] < tol) || (d1[k] > 1 - tol) ||
          (d2[k] < tol) || (d2[k] > 1 - tol) ||
          (1 - d1[k] - d2[k] < tol) || (1 - d1[k] - d2[k] > 1 - tol)) {

        /* lambda is set to ensure the lowest probability is MACHINE_TOL. */
        d1[k] = (1 - 3 * tol) * d1[k] + (3 * tol) * 1/3;
        d2[k] = (1 - 3 * tol) * d2[k] + (3 * tol) * 1/3;

      }/*THEN*/

    }/*THEN*/

    /* move to the next arc. */
    k++;

  }/*FOR*/

  /* set up the return value. */
  PROTECT(result = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(result, 0, from);
  SET_VECTOR_ELT(result, 1, to);
  SET_VECTOR_ELT(result, 2, nid);
  SET_VECTOR_ELT(result, 3, dir1);
  SET_VECTOR_ELT(result, 4, dir2);
  setAttrib(result, R_NamesSymbol,
    mkStringVec(5, "from", "to", "aid", "fwd", "bkwd"));
  PROTECT(df = minimal_data_frame(result));

  Free1D(poset);

  UNPROTECT(11);

  return df;

}/*CASTELO_COMPLETION*/

double marginal_prior(SEXP beta, SEXP target, SEXP parents, SEXP children,
    SEXP node_cache, SEXP nodes, bool debugging) {

int i = 0, t = 0, nnodes = length(nodes);
int *temp = NULL;
double prior = 0, result = 0, b = NUM(beta);
short int *adjacent = NULL;
SEXP try;

  /* match the target node. */
  PROTECT(try = match(nodes, target, 0));
  t = INT(try);
  UNPROTECT(1);

  /* find out which nodes are parents and which nodes are children. */
  adjacent = Calloc1D(nnodes, sizeof(short int));

  PROTECT(try = match(nodes, parents, 0));
  temp = INTEGER(try);
  for (i = 0; i < length(try); i++)
    adjacent[temp[i] - 1] = PARENT;
  UNPROTECT(1);

  PROTECT(try = match(nodes, children, 0));
  temp = INTEGER(try);
  for (i = 0; i < length(try); i++)
    adjacent[temp[i] - 1] = CHILD;
  UNPROTECT(1);

  /* prior probabilities table lookup. */
  for (i = t + 1; i <= nnodes; i++) {

    /* look up the prior probability. */
    prior = (adjacent[i - 1] > 0) ? b / 2 : 1 - b;

    if (debugging) {

      switch(adjacent[i - 1]) {

        case PARENT:
          Rprintf("  > found arc %s -> %s, prior pobability is %lf.\n",
            NODE(i - 1), NODE(t - 1), prior);
          break;
        case CHILD:
          Rprintf("  > found arc %s -> %s, prior probability is %lf.\n",
            NODE(t - 1), NODE(i - 1), prior);
          break;
        default:
          Rprintf("  > no arc between %s and %s, prior probability is %lf.\n",
            NODE(t - 1), NODE(i - 1), prior);

      }/*SWITCH*/

    }/*THEN*/

    /* move to log-scale and divide by the non-informative log(1/3), so that
     * the contribution of each arc whose prior has not been specified by the
     * user is zero; overflow is likely otherwise. */
    result += log(prior / ((double)1/3));

  }/*FOR*/

  Free1D(adjacent);

  return result;

}/*MARGINAL_PRIOR*/

