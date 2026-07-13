#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/data.table.h"
#include "../../fitted/fitted.h"
#include "../../include/globals.h"
#include "../../include/sampling.h"
#include "../../minimal/common.h"
#include "../../minimal/data.frame.h"
#include "../predict.h"

/* predict the values of one or more variables given one or more variables by
 * maximum a posteriori (MAP). */
SEXP predict_map_lw(SEXP node, SEXP fitted, SEXP data, SEXP n, SEXP from, SEXP prob,
    SEXP debug) {

int j = 0, nobs = 0, nev = 0, nlvls = 0, drop = 0;
int node_id = 0, nsims = INT(n);
bool *ev_discrete = NULL;
fixed_node *fixed = NULL;
void **varptrs = NULL, **evptrs = NULL, *pred = NULL, *res = NULL;
SEXP result, colnames, evidence, evmatch, temp = R_NilValue;
SEXP cpdist, predicted, lvls = R_NilValue, probtab = R_NilValue;
double *wgt = NULL, *pt = NULL;
long double *lvls_counts = NULL;
bool debugging = isTRUE(debug), include_prob = isTRUE(prob);
bool target_discrete = FALSE;
fitted_bn bn = fitted_network_from_SEXP(fitted);
tabular dt = tabular_from_SEXP(data, 0, 0);

  /* find the index of the target node in the fitted network, and whether it is
   * a discrete node. */
  for (node_id = 0; node_id < bn.nnodes; node_id++)
    if (strcmp(bn.labels[node_id], CHAR(STRING_ELT(node, 0))) == 0)
      break;
  target_discrete = (bn.node_types[node_id] == DNODE) ||
                      (bn.node_types[node_id] == ONODE);

  /* extract the names of the variables in the data. */
  colnames = getAttrib(data, R_NamesSymbol);

  /* remove the name of the variable to predict. */
  nev = length(from);
  PROTECT(evmatch = match(colnames, from, 0));

  /* cache, for each evidence variable, whether it is discrete and a pointer to
   * its column in the unified data table. */
  ev_discrete = Calloc1D(nev, sizeof(bool));
  varptrs = (void **) Calloc1D(nev, sizeof(void *));
  for (j = 0; j < nev; j++) {

    int dc = INTEGER(evmatch)[j] - 1;
    ev_discrete[j] = dt.m.flag[dc].discrete;
    varptrs[j] = ev_discrete[j] ?
                   (void *) dt.dcol[dt.map[dc]] : (void *) dt.ccol[dt.map[dc]];

  }/*FOR*/

  /* cache the sample size. */
  nobs = dt.m.nobs;

  /* allocate a list to hold the evidence. */
  PROTECT(evidence = allocVector(VECSXP, nev));
  setAttrib(evidence, R_NamesSymbol, from);

  /* cache pointers to the elements of the evidence .*/
  evptrs = (void **) Calloc1D(nev, sizeof(void *));

  for (j = 0; j < nev; j++) {

    PROTECT(temp = allocVector(ev_discrete[j] ? INTSXP : REALSXP, 1));
    evptrs[j] = DATAPTR(temp);
    SET_VECTOR_ELT(evidence, j, temp);
    UNPROTECT(1);

  }/*FOR*/

  /* convert the evidence and point evptrs into the the evidence. */
  fixed = evidence_from_SEXP(evidence, bn);
  for (j = 0; j < nev; j++) {

    int node = 0;

    for (node = 0; node < bn.nnodes; node++)
      if (strcmp(bn.labels[node], CHAR(STRING_ELT(from, j))) == 0)
        break;

    evptrs[j] = ev_discrete[j] ?
                  (void *) fixed[node].idx : (void *) fixed[node].val;

  }/*FOR*/

  /* allocate the return value. */
  PROTECT(result = fitnode2df(fitted, STRING_ELT(node, 0), nobs));
  res = DATAPTR(result);

  if (target_discrete) {

    /* for discrete variables, allocate scratch space for levels' frequencies
     * and for the prediction probabilities (if needed). */
    lvls = getAttrib(result, R_LevelsSymbol);
    nlvls = bn.ldists[node_id].d.dims[0];
    lvls_counts = Calloc1D(nlvls, sizeof(long double));

    if (include_prob) {

      PROTECT(probtab = allocMatrix(REALSXP, nlvls, nobs));
      pt = REAL(probtab);
      memset(pt, '\0', nobs * nlvls * sizeof(double));

    }/*THEN*/

  }/*THEN*/
  else {

    /* for continuous variables, there are no prediction probabilities. */
    include_prob = FALSE;

  }/*ELSE*/

  /* allocate the weights. */
  wgt = Calloc1D(nsims, sizeof(double));

  /* allocate sratch space for the random samplings. */
  PROTECT(cpdist = fit2df(fitted, nsims));
  predicted = getListElement(cpdist, (char *)CHAR(STRING_ELT(node, 0)));
  pred = DATAPTR(predicted);

  /* run the predictions on the bnlearn data structures. */
  c_mappred(dt, bn, node_id, target_discrete, nlvls, nev, ev_discrete, varptrs,
    evptrs, wgt, nsims, lvls_counts, pred, res, pt, include_prob, &drop,
    debugging, fitted, cpdist, from, fixed);

  /* deallocate here to avoid leaking memory if warnings are errors. */
  Free1D(ev_discrete);
  Free1D(varptrs);
  Free1D(evptrs);
  Free1D(wgt);
  if (target_discrete)
    Free1D(lvls_counts);
  FreeTAB(dt);
  FreeEvidence(fixed, bn.nnodes);
  FreeFittedBN(bn);

  if (drop > 0)
    warning("dropping %d observations because generated samples are NAs.", drop);

  if (include_prob) {

    /* set the levels of the taregt variable as rownames. */
    setDimNames(probtab, lvls, R_NilValue);
    /* add the posterior probabilities to the return value. */
    setAttrib(result, BN_ProbSymbol, probtab);

    UNPROTECT(5);

  }/*THEN*/
  else {

    UNPROTECT(4);

  }/*ELSE*/

  return result;

}/*PREDICT_MAP_LW*/
