#include "include/rcore.h"
#include "include/tests.h"

double discrete_df(test_e test, int *ni, int llx, int *nj, int lly) {

int i = 0, j = 0, alx = 0, aly = 0;
double df = 0;

  switch(test) {

    /* usual degrees of freedom. */
    case MI:
    case MI_SH:
    case X2:
      df = (llx - 1) * (lly - 1);
      break;

    /* adjust degrees of freedom: zeroes are considered structural if they
     * are part of a column or a row with a zero marginal, the rest are
     * considered sampling zeros. */
    case MI_ADF:
    case X2_ADF:
      for (i = 0; i < llx; i++)
        alx += (ni[i] > 0);
      for (j = 0; j < lly; j++)
        aly += (nj[j] > 0);

      /* ensure the degrees of freedom will not be negative. */
      alx = (alx >= 1) ? alx : 1;
      aly = (aly >= 1) ? aly : 1;

      df = (alx - 1) * (aly - 1);

      break;

    default:
      error("no degrees of freedom for this test.");

  }/*SWITCH*/

  return df;

}/*DISCRETE_DF*/

double discrete_cdf(test_e test, int **ni, int llx, int **nj, int lly, int llz) {

int i = 0, j = 0, k = 0, alx = 0, aly = 0;
double df = 0;

  switch(test) {

    /* usual degrees of freedom. */
    case MI:
    case MI_SH:
    case X2:
      df = (llx - 1) * (lly - 1) * llz;
      break;

    /* adjust degrees of freedom: zeroes are considered structural if they
     * are part of a column or a row with a zero marginal, the rest are
     * considered sampling zeros. */
    case MI_ADF:
    case X2_ADF:
      for (k = 0; k < llz; k++) {

        alx = aly = 0;

        for (i = 0; i < llx; i++)
          alx += (ni[k][i] > 0);
        for (j = 0; j < lly; j++)
          aly += (nj[k][j] > 0);

        /* ensure the degrees of freedom will not be negative. */
        alx = (alx >= 1) ? alx : 1;
        aly = (aly >= 1) ? aly : 1;

        df += (alx - 1) * (aly - 1);

      }/*FOR*/

      break;

    default:
      error("no degrees of freedom for this test.");

  }/*SWITCH*/

  return df;

}/*DISCRETE_CDF*/

double gaussian_cdf(test_e test, int num, int nz) {

double df = 0;

  switch(test) {

    case COR:
      df = num - nz - 2;
      break;

    case MI_G:
    case MI_G_SH:
      df = 1;
      break;

    case ZF:
      df = num - nz - 3;
      break;

    default:
      error("no degrees of freedom for this test.");

  }/*SWITCH*/

  return df;


}/*GAUSSIAN_CDF*/
