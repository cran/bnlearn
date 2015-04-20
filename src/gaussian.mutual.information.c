#include "include/rcore.h"
#include "include/covariance.h"

/* unconditional Gaussian mutual information, to be used in C code. */
double c_mig(double *xx, double *yy, int *num) {

double cor = c_fast_cor(xx, yy, *num);

  return - 0.5 * log(1 - cor * cor);

}/*C_MIG*/

