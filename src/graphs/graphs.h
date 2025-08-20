#ifndef GRAPHICAL_FUNCTIONS_HEADER
#define GRAPHICAL_FUNCTIONS_HEADER

int ug_connected_components(int *amat, char **labels, int nnodes, int **buffer,
    int *buflen, bool debugging);
void dag_path_matrix(int *amat, int *pathmat, int nnodes, char **labels,
    int *root_ids, int nroots, bool debugging);
void initialize_reachability_matrix(int *amat, int *pmat, int nnodes,
    int target, bool *cond_set, int *rmat, char **labels, bool debugging);
void complete_reachability_matrix(int *initial_rmat, int *final_rmat,
    int nnodes, char **labels, bool debugging);

#endif
