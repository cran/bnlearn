#include "../include/rcore.h"
#include "../core/covariance.matrix.h"
#include "../math/linear.algebra.h"

/* check that a matrix is symmetric and satisfies Cauchy-Schwarz. */
SEXP check_covariance(SEXP covmat) {

int i = 0, j = 0, n = sqrt(length(covmat));
double *cov = REAL(covmat);

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++) {

      /* firstly, check symmetry. */
      if (cov[CMC(i, j, n)] != cov[CMC(j, i, n)])
        error("'covmat' (%d, %d) is not symmetric.", i + 1, j + 1);

      /* secondly, check Cauchy-Schwarz. */
      if (cov[CMC(i, j, n)] > sqrt(cov[CMC(i, i, n)] * cov[CMC(j, j, n)]))
        error("'covmat' (%d, %d) does not satisfy the Cauchy-Schwarz "
          "inequality.", i + 1, j + 1);

    }/*FOR*/

  return R_NilValue;

}/*CHECK_COVARIANCE*/
