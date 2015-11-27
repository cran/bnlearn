#include "include/rcore.h"
#include "include/covariance.h"

/* correlation transform for the exact t-test. */
double cor_t_trans(double cor, double df) {

  return cor * sqrt(df) / sqrt(1 - cor * cor);

}/*COR_T_TRANS*/

/* correlation transform for Fizher's Z test. */
double cor_zf_trans(double cor, double df) {

  return log((1 + cor)/(1 - cor)) / 2 * sqrt(df - 1);

}/*COR_ZF_TRANS*/

/* correlation transform for mutual information. */
double cor_mi_trans(double cor) {

  return - 0.5 * log(1 - cor * cor);

}/*COR_MI_TRANS*/

