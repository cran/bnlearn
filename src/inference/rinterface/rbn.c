#include "../../include/rcore.h"
#include "../../core/allocations.h"
#include "../../core/data.table.h"
#include "../../fitted/fitted.h"
#include "../../include/sampling.h"
#include "../../minimal/data.frame.h"
#include "../../sanitization/data.h"

/* build the per-node evidence structure from the SEXP. */
fixed_node *evidence_from_SEXP(SEXP fix, fitted_bn fitted) {

fixed_node *fixed = Calloc1D(fitted.nnodes, sizeof(fixed_node));
SEXP fix_names;

  /* no evidence: fix is a logical (FALSE). */
  if (TYPEOF(fix) == LGLSXP)
    return fixed;

  fix_names = getAttrib(fix, R_NamesSymbol);

  for (int j = 0; j < length(fix); j++) {

    /* find the node this evidence entry refers to. */
    const char *name = CHAR(STRING_ELT(fix_names, j));
    int node = -1;

    for (int i = 0; i < fitted.nnodes; i++)
      if (strcmp(fitted.labels[i], name) == 0) {

        node = i;
        break;

      }/*THEN*/

    if (node < 0)
      continue;

    SEXP el = VECTOR_ELT(fix, j);
    int n = length(el);

    fixed[node].fixed = TRUE;
    fixed[node].n = n;

    if ((fitted.node_types[node] == DNODE) ||
        (fitted.node_types[node] == ONODE)) {

      /* discrete: store the 1-based level indexes. */
      fixed[node].idx = Calloc1D(n, sizeof(int));

      if (TYPEOF(el) == INTSXP) {

        for (int k = 0; k < n; k++)
          fixed[node].idx[k] = INTEGER(el)[k];

      }/*THEN*/
      else {

        char **levels = fitted.ldists[node].d.levels;
        int nlevels = fitted.ldists[node].d.dims[0];

        for (int k = 0; k < n; k++) {

          const char *label = CHAR(STRING_ELT(el, k));

          for (int l = 0; l < nlevels; l++)
            if (strcmp(levels[l], label) == 0) {

              fixed[node].idx[k] = l + 1;
              break;

            }/*THEN*/

        }/*FOR*/

      }/*ELSE*/

    }/*THEN*/
    else {

      /* continuous: store the value(s) (point or interval). */
      fixed[node].val = Calloc1D(n, sizeof(double));
      memcpy(fixed[node].val, REAL(el), n * sizeof(double));

    }/*ELSE*/

  }/*FOR*/

  return fixed;

}/*EVIDENCE_FROM_SEXP*/

void FreeEvidence(fixed_node *fixed, int nnodes) {

  for (int i = 0; i < nnodes; i++) {

    Free1D(fixed[i].idx);
    Free1D(fixed[i].val);

  }/*FOR*/

  Free1D(fixed);

}/*FREEEVIDENCE*/

/* generate random observations from a bayesian network. */
SEXP rbn_master(SEXP fitted, SEXP n, SEXP fix, SEXP add_metadata, SEXP debug) {

int num = INT(n);
fitted_bn bn = fitted_network_from_SEXP(fitted);
fixed_node *fixed = evidence_from_SEXP(fix, bn);
tabular result_tab = { 0 };
SEXP result;

  /* allocate the return value. */
  PROTECT(result = fit2df(fitted, num));
  /* build a tabular view of the result and generate the random observations. */
  result_tab = tabular_from_SEXP(result, 0, 0);
  c_rbn_master(bn, result_tab, fixed, isTRUE(debug));

  /* add the metadata, if requested. */
  if (isTRUE(add_metadata))
    add_simulation_metadata(result, num);

  FreeTAB(result_tab);
  FreeEvidence(fixed, bn.nnodes);
  FreeFittedBN(bn);

  UNPROTECT(1);

  return result;

}/*RBN_MASTER*/
