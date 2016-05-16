
/* from filter.arcs.c */
SEXP which_undirected(SEXP arcs, SEXP nodes);

/* from tiers.c */
SEXP tiers(SEXP nodes, SEXP debug);

/* from filter.arcs.c */
SEXP c_unique_arcs(SEXP arcs, SEXP nodes, int warnlevel);

/* from fitted.c */
SEXP root_nodes(SEXP bn, SEXP leaves);

/* from simulation.c */
SEXP schedule(SEXP bn, SEXP root_nodes, SEXP reverse, SEXP debug);

/* from path.c */
int c_has_path(int start, int stop, int *amat, int n, SEXP nodes,
    int ugraph, int notdirect,  int *path, int *counter, int debuglevel);
int c_directed_path(int start, int stop, int *amat, int n, SEXP nodes,
    int *path, int *counter, int debuglevel);
int c_uptri3_path(short int *uptri, int *depth, int from, int to, int nnodes,
    SEXP nodes, int debuglevel);

/* from hash.c */
SEXP arc_hash(SEXP arcs, SEXP nodes, int uptri, int sort);
void c_arc_hash(int narcs, int nnodes, int *from, int *to, int *uptri,
  int *cmc, int sort);
SEXP c_amat_hash(int *amat, int nnodes);

/* from arcs2amat.c */
SEXP arcs2amat(SEXP arcs, SEXP nodes);
SEXP amat2arcs(SEXP amat, SEXP nodes);

