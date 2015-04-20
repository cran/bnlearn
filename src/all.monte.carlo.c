#include "include/rcore.h"
#include "include/tests.h"

int remap_permutation_test(const char *t) {

  /* remap the test statistics to the constants used in monte.carlo.c. */
  if ((strcmp(t, "mc-mi") == 0) || (strcmp(t, "smc-mi") == 0))
    return MUTUAL_INFORMATION;
  else if ((strcmp(t, "mc-x2") == 0) || (strcmp(t, "smc-x2") == 0))
    return PEARSON_X2;
  else if ((strcmp(t, "mc-mi-g") == 0) || (strcmp(t, "smc-mi-g") == 0))
    return GAUSSIAN_MUTUAL_INFORMATION;
  else if ((strcmp(t, "mc-cor") == 0) || (strcmp(t, "smc-cor") == 0))
    return LINEAR_CORRELATION;
  else if ((strcmp(t, "mc-zf") == 0) || (strcmp(t, "smc-zf") == 0))
    return FISHER_Z;
  else if (strcmp(t, "sp-mi") == 0)
    return SP_MUTUAL_INFORMATION;
  else if (strcmp(t, "sp-x2") == 0)
    return SP_PEARSON_X2;
  else if ((strcmp(t, "mc-jt") == 0) || (strcmp(t, "smc-jt") == 0))
    return JT;

  return -1;

}/*REMAP_PERMUTATION_TEST*/

