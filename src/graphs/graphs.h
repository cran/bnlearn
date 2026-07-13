#ifndef GRAPHICAL_FUNCTIONS_HEADER
#define GRAPHICAL_FUNCTIONS_HEADER

#include "../math/sparse.amat.h"
#include "../core/data.table.h"
#include "../parameters/parameters.h"

int ug_connected_components(int *amat, char **labels, int nnodes, int **buffer,
    int *buflen, bool debugging);

/* number of parameters of a single node (from nparams.c). */
double nparams_node_count(sparse_amat parents, tabular dt, int i,
    estimator_e estimator);

#endif
