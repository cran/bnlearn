#include "../core/contingency.tables.h"

/* from rbn.c */
void c_rbn_master(SEXP fitted, SEXP result, SEXP n, SEXP fix, bool add_metadata,
    bool debugging);

/* from likelihood.weighting.c */
void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    bool debugging);
