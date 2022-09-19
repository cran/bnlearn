#include "../../include/rcore.h"
#include "../tests.h"

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
