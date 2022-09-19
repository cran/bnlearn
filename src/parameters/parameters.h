#ifndef PARAMETER_LEARNING_HEADER
#define PARAMETER_LEARNING_HEADER

#include "../math/linear.algebra.h"

/* bit-field tracking all possible errors in the hierarchical Dirichlet
 * parameter estimation. */
typedef struct {

  unsigned int outer_em_convergence_fail  : 1;
  unsigned int kappa_tau_convergence_fail : 1;
  unsigned int tau_convergence_fail       : 1;
  unsigned int tau_is_zero                : 1;
  unsigned int kappa_convergence_fail     : 1;
  unsigned int padding                    : 3;  /* pad to 1 byte. */

} hdstatus;

void c_classic_discrete_parameters(int *counts, double *cpt, int nrows,
    int ncols, double alpha, bool replace);
hdstatus c_hierarchical_dirichlet_parameters(cmcmap counts, double alpha0,
    double s, bool debugging, double *nu);

#endif

