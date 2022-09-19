#include "../include/rcore.h"
#include "../minimal/common.h"

/* Kullback-Leibler divergence between conditional probability tables. */
SEXP kullback_leibler_discrete(SEXP cptableP, SEXP cptableQ) {

double *p = REAL(cptableP), *q = REAL(cptableQ);
long double divergence = 0;

  /* the convention is 0 log 0 = 0 and 0 log 0 / 0 = 0. */
  for (int i = 0; i < length(cptableP); i++)
    if (p[i] != 0)
      divergence += p[i] * log(p[i] / q[i]);

  return mkReal(divergence);

}/*KULLBACK_LEIBLER_DISCRETE*/

