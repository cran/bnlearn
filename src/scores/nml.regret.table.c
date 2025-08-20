#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../math/linear.algebra.h"
#include "scores.h"

double get_regret_k2(int N, double *logs, double *logfacs) {

double res = 0.0;

  for (int n = 0; n <= N; ++n) {

    double term = logfacs[N] - logfacs[n] - logfacs[N - n];

    if (n)
      term += n * (logs[n] - logs[N]);
    if (N - n)
      term += (N - n) * (logs[N - n] - logs[N]);

    res += exp(term);

  }/*FOR*/

  return log(res);

}/*GET_REGRET_K2*/

void fill_regrets_up_to_K(int K, int n, double *regret_table) {

double rk_prev_prev = 1.0;
double rk_prev = exp(regret_table[CMC(2, n, K + 1)]);
double rk = 0;

  for (int k = 3; k <= K; ++k) {

    rk = rk_prev + rk_prev_prev / (k - 2) * n;
    regret_table[CMC(k, n, K + 1)] = log(rk);
    rk_prev_prev = rk_prev;
    rk_prev = rk;

  }/*FOR*/

}/*FILL_REGRETS_UP_TO_K*/

double *get_regret_table(int N, int K) {

double *logs = (double *) Calloc1D((N + 1), sizeof(double));
double *logfacs = (double *) Calloc1D((N + 1), sizeof(double));
double *regret_table = (double *) Calloc1D((N + 1) * (K + 1), sizeof(double));

  /* cache a table of logarithms. */
  for (int n = 1; n <= N; ++n)
    logs[n] = log(n);

  /* cache a table of log-factorials (computed as log-Gammas). */
  for (int n = 1; n <= N; ++n)
    logfacs[n] = lgammafn(n + 1);

  for (int n = 1; n <= N; ++n) {

    regret_table[CMC(2, n, K + 1)] = get_regret_k2(n, logs, logfacs);
    fill_regrets_up_to_K(K, n, regret_table);

  }/*FOR*/

  Free1D(logfacs);
  Free1D(logs);

  return regret_table;

}/*GET_REGRET_TABLE*/
