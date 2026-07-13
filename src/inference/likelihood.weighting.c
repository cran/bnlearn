#include "../include/rcore.h"
#include "../core/data.table.h"
#include "../core/math.functions.h"
#include "../fitted/fitted.h"
#include "loglikelihood/loglikelihood.h"

void c_lw_weights(fitted_bn bn, tabular dt, double *w, bool debugging) {

int i = 0, max_el = 0, n = dt.m.nobs;
double maxw = 0;

  memset(w, '\0', n * sizeof(double));

  /* compute the log-likelihood of each sample under the fitted network. */
  if (bn.type == DNET || bn.type == ONET || bn.type == DONET)
    bysample_discrete_loglikelihood(bn, dt, w, debugging);
  else if (bn.type == GNET)
    bysample_gaussian_loglikelihood(bn, dt, w, TRUE, debugging);
  else if (bn.type == CGNET)
    bysample_clgaussian_loglikelihood(bn, dt, w, TRUE, debugging);
  else if (bn.type == ZINET)
    bysample_zeroinflated_loglikelihood(bn, dt, w, TRUE, debugging);

  /* rescale before exponentiating into probabilities (if possible). */
  max_el = d_which_max(w, n);

  if (max_el == NA_INTEGER)
    return;
  else if ((max_el == 1) && (w[0] == R_NegInf))
    memset(w, '\0', n * sizeof(double));
  else {

    maxw = w[max_el - 1];
    for (i = 0; i < n; i++)
      w[i] = exp(w[i] - maxw);

  }/*ELSE*/

}/*C_LW_WEIGHTS*/
