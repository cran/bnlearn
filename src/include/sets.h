
/* from configurations.c */
void cfg(SEXP parents, int *configurations, int *nlevels);
void c_fast_config(int **columns, int nrows, int ncols, int *levels,
    int *configurations, int *nlevels, int offset);
SEXP c_configurations(SEXP parents, int factor, int all_levels);

/* from subsets.c */
void first_subset(int *work, int n, int offset);
int next_subset(int *work, int n, int max, int offset);

