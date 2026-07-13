#ifndef INFERENCE_SAMPLING_HEADER
#define INFERENCE_SAMPLING_HEADER

#include "../core/contingency.tables.h"
#include "../fitted/fitted.h"
#include "../core/data.table.h"

/* evidence (fixed nodes) for random generation, indexed by node. */
typedef struct {

  bool fixed;     /* whether this node has evidence. */
  int n;          /* number of values. */
  int *idx;       /* discrete: 1-based level indexes (length n), NULL otherwise. */
  double *val;    /* continuous: values (length n, 1 = point, 2 = interval),
                   * NULL otherwise. */

} fixed_node;

/* from rbn.c */
void c_rbn_master(fitted_bn fitted, tabular result, fixed_node *fixed,
    bool debugging);

/* from rinterface/rbn.c */
fixed_node *evidence_from_SEXP(SEXP fix, fitted_bn fitted);
void FreeEvidence(fixed_node *fixed, int nnodes);

/* from likelihood.weighting.c */
void c_lw_weights(fitted_bn fitted, tabular data, double *w, bool debugging);

/* from rinterface/likelihood.weighting.c */
tabular lw_data_from_SEXP(SEXP data, SEXP fitted, SEXP keep);

#endif
