
/* from cache.structure.c */
SEXP cache_structure(SEXP nodes, SEXP amat, SEXP debug);
SEXP c_cache_partial_structure(int target, SEXP nodes, SEXP amat, int *status,
    SEXP debug);

/* from common.c */
SEXP finalize_arcs(SEXP arcs);
