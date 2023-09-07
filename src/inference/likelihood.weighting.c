#include "../include/rcore.h"
#include "../core/math.functions.h"
#include "../inference/loss.h"

void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    bool debugging) {

int i = 0, max_el = 0;
double maxw = 0;

  /* ensure the buffer is clean. */
  memset(w, '\0', n * sizeof(double));
  /* compute log-probabilities for each particle. */
  c_entropy_loss(fitted, data, n, TRUE, w, NULL, keep, FALSE, FALSE, debugging);
  /* rescale before exponentiating them into probabilities (if possible). */
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

