#include "include/rcore.h"
#include "include/matrix.h"
#include "include/scores.h"

/* regret approximation for large N */
double regret_fn_szp1(double N, double K) {

double hK   = 0.5 * K;
double hK_h = hK - 0.5;
double lgammahK = lgamma(hK);
double lgammahK_h = lgamma(hK_h);
double gamma_ratio = exp(lgammahK - lgammahK_h);
double gamma_ratio_term = M_SQRT1_2 - gamma_ratio * K/3/sqrt(N);

  return (hK_h) * (log(N) - M_LN2)
           + M_LN_SQRT_PI - lgammahK
           + 0.5 - gamma_ratio_term * gamma_ratio_term
           + (3 + K * (K - 2) * (2 * K + 1)) / (36.0 * N);

}/*REGRET_FN_SZP1*/

/* regret approximation for large K */
double regret_fn_szp2(double N, double K) {

double a = K/N;
double Ca = 0.5 * (1 + sqrt(1 + 4/a));

  return N * (log(a) + (a + 2) * log(Ca) - 1/Ca) - 0.5 * log(Ca + 2/a);

}/*REGRET_FN_SZP2*/

/* allocate and fill the regret table. */
double nml_regret(double n, double k) {

  if (n == 1) {

    return log(k);

  }/*THEN*/
  else if (k == 1 || n == 0) {

    return 0.0;

  }/*THEN*/
  else if (n <= MAX_REGRET_TABLE_N && k <= MAX_REGRET_TABLE_K) {

    if (regret_table == NULL)
      regret_table = get_regret_table(MAX_REGRET_TABLE_N, MAX_REGRET_TABLE_K);

    return regret_table[CMC((int) k, (int) n, MAX_REGRET_TABLE_K + 1)];

  }/*THEN*/
  else {

    double approx1 = regret_fn_szp1(n, k);
    double approx2 = regret_fn_szp2(n, k);

    return (approx1 < approx2) ? approx1 : approx2;

  }/*ELSE*/

}/*NML_REGRET*/
