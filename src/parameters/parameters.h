#ifndef PARAMETER_LEARNING_HEADER
#define PARAMETER_LEARNING_HEADER

#include "../math/linear.algebra.h"

/* enum for the parameter estimators, to be matched from the label string passed
 * down from R (the available.fits vector in R/globals.R). */
typedef enum {
  ENOEST      =  0, /* error code, no such estimator. */

  /* estimators for discrete data. */
  MLE         =  1, /* maximum likelihood. */
  BAYES       =  2, /* Bayesian (posterior) Dirichlet. */
  HDIR        =  3, /* hierarchical Dirichlet. */
  HARD_EM     =  4, /* hard expectation-maximization. */

  /* estimators for Gaussian data. */
  MLE_G       =  5, /* maximum likelihood. */
  HARD_EM_G   =  6, /* hard expectation-maximization. */

  /* estimators for conditional Gaussian (mixed) data. */
  MLE_CG      =  7, /* maximum likelihood. */
  HARD_EM_CG  =  8, /* hard expectation-maximization. */

  /* estimators for zero-inflated count data. */
  MLE_ZIHP    =  9, /* zero-inflated hyper-Poisson. */
  MLE_ZINB    = 10  /* zero-inflated negative binomial. */
} estimator_e;

estimator_e estimator_to_enum(const char *label);

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

