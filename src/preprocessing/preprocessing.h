#ifndef PREPROCESSING_HEADER
#define PREPROCESSING_HEADER

#include "../core/data.table.h"

/* enum for discretization methods, to be matched from the label string passed
 * down from R. */
typedef enum {
  ENOMETHOD    =  0, /* error code, no such discretization method. */
  INTERVAL     =  1, /* interval discretization (marginal). */
  QUANTILE     =  2, /* quantile discretization (marginal). */
  HARTEMINK    =  3  /* Hartemink discretization (pairwise). */
} discretization_e;

discretization_e discretization_to_enum(const char *label);

int interval_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs, bool debugging);
int quantile_discretization(double *orig, int *factor, int nbreaks,
    double *cutpoints, int nobs, bool complete, bool debugging);
void hartemink_discretization(ddata work, int *nbreaks, double **cutpoints,
    bool debugging);

#endif
