#include "common.h"

/* unconditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int i = 0, ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
void *memp = NULL;
SEXP xx, yy, result;

  /* allocate the return value, which has the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  nobs = length(yy);

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    int llx = 0, lly = NLEVELS(yy), adj = (strcmp(t, "mi-adf") == 0);
    int *xptr = NULL, *yptr = INTEGER(yy);
    SEXP xdata;

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = 2 * nobs * c_mi(xptr, llx, yptr, lly, nobs, &df, adj);

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    int llx = 0, lly = NLEVELS(yy);
    int *xptr = NULL, *yptr = INTEGER(yy);
    SEXP xdata;

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = 2 * nobs * c_shmi(xptr, llx, yptr, lly, nobs);
      df = ((double)(llx - 1) * (double)(lly - 1));

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    int llx = 0, lly = NLEVELS(yy), adj = (strcmp(t, "x2-adf") == 0);
    int *xptr = NULL, *yptr = INTEGER(yy);
    SEXP xdata;

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = c_x2(xptr, llx, yptr, lly, nobs, &df, adj);

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    int llx = 0, lly = NLEVELS(yy);
    int *xptr = NULL, *yptr = INTEGER(yy);
    SEXP xdata;

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = c_jt(xptr, llx, yptr, lly, nobs);

      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* linear correlation. */
    double transform = 0, *xptr = NULL, *yptr = REAL(yy);
    double xm = 0, ym = 0, xsd = 0, ysd = 0;
    df = nobs - 2;

    /* check that we have enough degrees of freedom. */
    if (df < 1)
      error("trying to do an independence test with zero degrees of freedom.");

    /* cache mean and variance. */
    ym = c_mean(yptr, nobs);
    ysd = c_var(yptr, ym, nobs);

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */
      xptr = REAL(VECTOR_ELT(xx, i));
      xm = c_mean(xptr, nobs);
      xsd = c_var(xptr, xm, nobs);
      statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
      transform = fabs(statistic * sqrt(df) / sqrt(1 - statistic * statistic));
      pvalue[i] = 2 * pt(transform, df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    double *xptr = NULL, *yptr = REAL(yy);
    double xm = 0, ym = 0, xsd = 0, ysd = 0;

    /* check that we have enough degrees of freedom. */
    if (nobs - 3 < 1)
      error("trying to do an independence test with zero degrees of freedom.");

    /* cache mean and variance. */
    ym = c_mean(yptr, nobs);
    ysd = c_var(yptr, ym, nobs);

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */
      xptr = REAL(VECTOR_ELT(xx, i));
      xm = c_mean(xptr, nobs);
      xsd = c_var(xptr, xm, nobs);
      statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
      statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                    sqrt((double)nobs - 3);
      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    double *xptr = NULL, *yptr = REAL(yy);
    double xm = 0, ym = 0, xsd = 0, ysd = 0;
    df = 1;

    /* cache mean and variance. */
    ym = c_mean(yptr, nobs);
    ysd = c_var(yptr, ym, nobs);

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */
      xptr = REAL(VECTOR_ELT(xx, i));
      xm = c_mean(xptr, nobs);
      xsd = c_var(xptr, xm, nobs);
      statistic = c_fast_cor2(xptr, yptr, nobs, xm, ym, xsd, ysd);
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    double *yptr = REAL(yy);
    df = 1;

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */
      statistic = c_fast_shcor(REAL(VECTOR_ELT(xx, i)), yptr, &nobs);
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0;
    double a = (strncmp(t, "smc-", 4) == 0) ? NUM(alpha) : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type)) {

      int *xptr = NULL, *yptr = INTEGER(yy);
      int llx = 0, lly = NLEVELS(yy);
      SEXP xdata;

      for (i = 0; i < ntests; i++) {

        memp = vmaxget();

        xdata = VECTOR_ELT(xx, i);
        xptr = INTEGER(xdata);
        llx = NLEVELS(xdata);

        statistic = 0;

        c_mcarlo(xptr, llx, yptr, lly, nobs, INT(B), &statistic, pvalue + i,
          a, type, &df);

        vmaxset(memp);

      }/*FOR*/

    }/*THEN*/
    else {

      /* continuous test. */
      double *yptr = REAL(yy);

      for (i = 0; i < ntests; i++) {

        memp = vmaxget();

        statistic = 0;

        c_gauss_mcarlo(REAL(VECTOR_ELT(xx, i)), yptr, nobs, INT(B), pvalue + i,
          a, type, &statistic);

        vmaxset(memp);

      }/*FOR*/

    }/*ELSE*/

  }/*THEN*/

  UNPROTECT(3);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*UTEST*/

/* conditional independence tests. */
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning) {

int i = 0, ntests = length(x), nobs = 0;
double *pvalue = NULL, statistic = 0, df = NA_REAL;
const char *t = CHAR(STRING_ELT(test, 0));
void *memp = NULL;
SEXP xx, yy, zz, result;

  /* allocate the return value, which the same length as x. */
  PROTECT(result = allocVector(REALSXP, ntests));
  setAttrib(result, R_NamesSymbol, x);
  pvalue = REAL(result);
  /* set all elements to zero. */
  memset(pvalue, '\0', ntests * sizeof(double));

  /* extract the variables from the data. */
  PROTECT(xx = c_dataframe_column(data, x, FALSE, FALSE));
  PROTECT(yy = c_dataframe_column(data, y, TRUE, FALSE));
  PROTECT(zz = c_dataframe_column(data, sx, FALSE, FALSE));
  nobs = length(yy);

  if ((strcmp(t, "mi") == 0) || (strcmp(t, "mi-adf") == 0)) {

    /* mutual information test, with and without df adjustments. */
    int llx = 0, lly = NLEVELS(yy), llz = 0, adj = (strcmp(t, "mi-adf") == 0);
    int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
    SEXP xdata, config;

    PROTECT(config = c_configurations(zz, TRUE, TRUE));
    zptr = INTEGER(config);
    llz = NLEVELS(config);

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = 2 * nobs * c_cmi(xptr, llx, yptr, lly, zptr, llz, 
                               nobs, &df, adj);

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else if (strcmp(t, "mi-sh") == 0) {

    /* shrinkage mutual information test. */
    int llx = 0, lly = NLEVELS(yy), llz = 0;
    int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
    SEXP xdata, config;

    PROTECT(config = c_configurations(zz, TRUE, TRUE));
    zptr = INTEGER(config);
    llz = NLEVELS(config);

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = 2 * nobs * c_shcmi(xptr, llx, yptr, lly, zptr, llz,
                               nobs, &df);

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else if ((strcmp(t, "x2") == 0) || (strcmp(t, "x2-adf") == 0)) {

    /* Pearson's X^2 test, with and without df adjustments. */
    int llx = 0, lly = NLEVELS(yy), llz = 0, adj = (strcmp(t, "x2-adf") == 0);
    int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
    SEXP xdata, config;

    PROTECT(config = c_configurations(zz, TRUE, TRUE));
    zptr = INTEGER(config);
    llz = NLEVELS(config);

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = c_cx2(xptr, llx, yptr, lly, zptr, llz, nobs, &df, adj);

      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else if (strcmp(t, "jt") == 0) {

    /* Jonckheere-Terpstra test. */
    int llx = 0, lly = NLEVELS(yy), llz = 0;
    int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
    SEXP xdata, config;

    PROTECT(config = c_configurations(zz, TRUE, TRUE));
    zptr = INTEGER(config);
    llz = NLEVELS(config);

    for (i = 0; i < ntests; i++) {

      memp = vmaxget();

      xdata = VECTOR_ELT(xx, i);
      xptr = INTEGER(xdata);
      llx = NLEVELS(xdata);

      statistic = c_cjt(xptr, llx, yptr, lly, zptr, llz, nobs);

      pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

      vmaxset(memp);

    }/*FOR*/

    UNPROTECT(1);

  }/*THEN*/
  else if (strcmp(t, "cor") == 0) {

    /* Pearson's linear correlation test. */
    int nsx = length(sx), ncols = nsx + 2;
    double transform = 0, **column = NULL, *mean = NULL;
    double *u = NULL, *d = NULL, *vt = NULL, *cov = NULL, *basecov = 0;
    df = nobs - ncols;

    /* check that we have enough degrees of freedom. */
    if (df < 1)
      error("trying to do a conditional independence test with zero degrees of freedom.");

    /* allocate the covariance matrix. */
    cov = alloc1dreal(ncols * ncols);
    basecov = alloc1dreal(ncols * ncols);

    /* allocate the matrices needed for the SVD decomposition. */
    u = alloc1dreal(ncols * ncols);
    d = alloc1dreal(ncols);
    vt = alloc1dreal(ncols * ncols);

    /* allocate and initialize an array of pointers for the variables. */
    column = (double **) alloc1dpointer(ncols);
    column[1] = REAL(yy);
    for (i = 0; i < nsx; i++)
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

    if (ntests > 1) {

      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 1);
      c_covmat(column, mean, ncols, nobs, cov, 1);
      memcpy(basecov, cov, ncols * ncols * sizeof(double));

      for (i = 0; i < ntests; i++) {

        /* no allocations require a vmax{get,set}() call. */

        /* extract and de-reference the i-th variable. */
        column[0] = REAL(VECTOR_ELT(xx, i));
        /* update the corresponding mean in the cache. */
        c_update_meanvec(column, mean, 0, nobs);
        /* update the covariance matrix. */
        memcpy(cov, basecov, ncols * ncols * sizeof(double));
        c_update_covmat(column, mean, 0, ncols, nobs, cov);
        /* compute the partial correlation and the test statistic. */
        statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
        transform = fabs(statistic * sqrt(df) / sqrt(1 - statistic * statistic));
        pvalue[i] = 2 * pt(transform, df, FALSE, FALSE);

      }/*FOR*/

    }/*THEN*/
    else {

      /* extract and de-reference the i-th variable. */
      column[0] = REAL(VECTOR_ELT(xx, 0));
      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 0);
      c_covmat(column, mean, ncols, nobs, cov, 0);
      /* compute the partial correlation and the test statistic. */
      statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
      transform = fabs(statistic * sqrt(df) / sqrt(1 - statistic * statistic));
      pvalue[0] = 2 * pt(transform, df, FALSE, FALSE);

    }/*ELSE*/

  }/*THEN*/
  else if (strcmp(t, "zf") == 0) {

    /* Fisher's Z test. */
    int nsx = length(sx), ncols = nsx + 2;
    double **column = NULL, *mean = NULL, *cov = NULL, *basecov = NULL;
    double *u = NULL, *d = NULL, *vt = NULL;

    /* check that we have enough degrees of freedom. */
    if (nobs - 1 - ncols < 1)
      error("trying to do an independence test with zero degrees of freedom.");

    /* allocate the covariance matrix. */
    cov = alloc1dreal(ncols * ncols);
    basecov = alloc1dreal(ncols * ncols);

    /* allocate the matrices needed for the SVD decomposition. */
    u = alloc1dreal(ncols * ncols);
    d = alloc1dreal(ncols);
    vt = alloc1dreal(ncols * ncols);

    /* allocate and initialize an array of pointers for the variables. */
    column = (double **) alloc1dpointer(ncols);
    column[1] = REAL(yy);
    for (i = 0; i < nsx; i++)
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

    if (ntests > 1) {

      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 1);
      c_covmat(column, mean, ncols, nobs, cov, 1);
      memcpy(basecov, cov, ncols * ncols * sizeof(double));

      for (i = 0; i < ntests; i++) {

        /* no allocations require a vmax{get,set}() call. */

        /* extract and de-reference the i-th variable. */
        column[0] = REAL(VECTOR_ELT(xx, i));
        /* update the corresponding mean in the cache. */
        c_update_meanvec(column, mean, 0, nobs);
        /* update the covariance matrix. */
        memcpy(cov, basecov, ncols * ncols * sizeof(double));
        c_update_covmat(column, mean, 0, ncols, nobs, cov);
        /* compute the partial correlation and the test statistic. */
        statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
        statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                      sqrt((double)nobs - 1 - ncols);
        pvalue[i] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

      }/*FOR*/

    }/*THEN*/
    else {

      /* extract and de-reference the i-th variable. */
      column[0] = REAL(VECTOR_ELT(xx, 0));
      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 0);
      c_covmat(column, mean, ncols, nobs, cov, 0);
      /* compute the partial correlation and the test statistic. */
      statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
      statistic = log((1 + statistic)/(1 - statistic)) / 2 *
                    sqrt((double)nobs - 1 - ncols);
      pvalue[0] = 2 * pnorm(fabs(statistic), 0, 1, FALSE, FALSE);

    }/*ELSE*/

  }/*THEN*/
  else if (strcmp(t, "mi-g") == 0) {

    /* Gaussian mutual information test. */
    int nsx = length(sx), ncols = nsx + 2;
    double **column = NULL, *mean = NULL, *cov = NULL, *basecov = NULL;
    double *u = NULL, *d = NULL, *vt = NULL;
    df = 1;

    /* allocate the covariance matrix. */
    cov = alloc1dreal(ncols * ncols);
    basecov = alloc1dreal(ncols * ncols);

    /* allocate the matrices needed for the SVD decomposition. */
    u = alloc1dreal(ncols * ncols);
    d = alloc1dreal(ncols);
    vt = alloc1dreal(ncols * ncols);

    /* allocate and initialize an array of pointers for the variables. */
    column = (double **) alloc1dpointer(ncols);
    column[1] = REAL(yy);
    for (i = 0; i < nsx; i++)
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

    if (ntests > 1) {

      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 1);
      c_covmat(column, mean, ncols, nobs, cov, 1);
      memcpy(basecov, cov, ncols * ncols * sizeof(double));

      for (i = 0; i < ntests; i++) {

        /* no allocations require a vmax{get,set}() call. */

        /* extract and de-reference the i-th variable. */
        column[0] = REAL(VECTOR_ELT(xx, i));
        /* update the corresponding mean in the cache. */
        c_update_meanvec(column, mean, 0, nobs);
        /* update the covariance matrix. */
        memcpy(cov, basecov, ncols * ncols * sizeof(double));
        c_update_covmat(column, mean, 0, ncols, nobs, cov);
        /* compute the partial correlation and the test statistic. */
        statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
        statistic = - nobs * log(1 - statistic * statistic);
        pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

      }/*FOR*/

    }/*THEN*/
    else {

      /* extract and de-reference the i-th variable. */
      column[0] = REAL(VECTOR_ELT(xx, 0));
      /* allocate and compute mean values and the covariance matrix. */
      mean = alloc1dreal(ncols);
      c_meanvec(column, mean, nobs, ncols, 0);
      c_covmat(column, mean, ncols, nobs, cov, 0);
      /* compute the partial correlation and the test statistic. */
      statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[0] = pchisq(statistic, df, FALSE, FALSE);

    }/*ELSE*/

  }/*THEN*/
  else if (strcmp(t, "mi-g-sh") == 0) {

    /* shrinkage Gaussian mutual information test. */
    int nsx = length(sx), ncols = nsx + 2;
    double **column = NULL, *mean = NULL, *cov = NULL;
    double *u = NULL, *d = NULL, *vt = NULL;
    df = 1;

    /* allocate the covariance matrix. */
    cov = alloc1dreal(ncols * ncols);

    /* allocate the matrices needed for the SVD decomposition. */
    u = alloc1dreal(ncols * ncols);
    d = alloc1dreal(ncols);
    vt = alloc1dreal(ncols * ncols);

    /* allocate and initialize an array of pointers for the variables. */
    column = (double **) alloc1dpointer(ncols);
    column[1] = REAL(yy);
    for (i = 0; i < nsx; i++)
      column[i + 2] = REAL(VECTOR_ELT(zz, i));

    /* allocate and compute the mean values. */
    mean = alloc1dreal(ncols);
    c_meanvec(column, mean, nobs, ncols, 1);

    for (i = 0; i < ntests; i++) {

      /* no allocations require a vmax{get,set}() call. */

      /* extract and de-reference the i-th variable. */
      column[0] = REAL(VECTOR_ELT(xx, i));
      /* update the corresponding mean in the cache. */
      c_update_meanvec(column, mean, 0, nobs);
      /* compute the covariance matrix. */
      memset(cov, '\0', ncols * ncols * sizeof(double));
      c_cov_lambda(column, mean, ncols, nobs, cov);
      /* compute the partial correlation and the test statistic. */
      statistic = c_fast_pcor(cov, u, d, vt, ncols, TRUE);
      statistic = - nobs * log(1 - statistic * statistic);
      pvalue[i] = pchisq(statistic, df, FALSE, FALSE);

    }/*FOR*/

  }/*THEN*/
  else if ((strncmp(t, "mc-", 3) == 0) || (strncmp(t, "smc-", 4) == 0) ||
           (strncmp(t, "sp-", 3) == 0)) {

    /* nonparametric and semiparametric tests. */
    int type = 0;
    double a = (strncmp(t, "smc-", 4) == 0) ? NUM(alpha) : 1;

    /* remap the test statistics to the constants used in monte.carlo.c. */
    type = remap_permutation_test(t);

    if (DISCRETE_PERMUTATION_TEST(type)) {

      int *xptr = NULL, *yptr = INTEGER(yy), *zptr = NULL;
      int llx = 0, lly = NLEVELS(yy), llz = 0;
      SEXP xdata, config;

      PROTECT(config = c_configurations(zz, TRUE, TRUE));
      zptr = INTEGER(config);
      llz = NLEVELS(config);

      for (i = 0; i < ntests; i++) {

        memp = vmaxget();

        xdata = VECTOR_ELT(xx, i);
        xptr = INTEGER(xdata);
        llx = NLEVELS(xdata);

        statistic = 0;

        c_cmcarlo(xptr, llx, yptr, lly, zptr, llz, nobs, INT(B), &statistic,
          pvalue + i, a, type, &df);

        vmaxset(memp);

      }/*FOR*/

      UNPROTECT(1);

    }/*THEN*/
    else {

      int nsx = length(sx), ncols = nsx + 2;
      double **column = NULL, *yptr = REAL(yy);

      /* allocate and initialize an array of pointers for the variables. */
      column = (double **) alloc1dpointer(ncols);
      column[1] = REAL(yy);
      for (i = 0; i < nsx; i++)
        column[i + 2] = REAL(VECTOR_ELT(zz, i));

      for (i = 0; i < ntests; i++) {

        memp = vmaxget();

        /* swap the first column and restore the second, which is that undergoing
         * permutation (backward compatibility from set random seed). */
        column[0] = REAL(VECTOR_ELT(xx, i));
        column[1] = yptr;

        statistic = 0;

        c_gauss_cmcarlo(column, ncols, nobs, INT(B), &statistic, pvalue + i,
          a, type);

        vmaxset(memp);

      }/*FOR*/

    }/*ELSE*/

  }/*THEN*/

  UNPROTECT(4);

  /* increase the test counter. */
  test_counter += ntests;

  if (isTRUE(learning))
    return result;
  else
    return c_create_htest(statistic, test, pvalue[ntests - 1], df, B);

}/*CTEST*/

/* independence tests, frontend to be used in R code. */
SEXP indep_test(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP learning) {

  /* if either node to test is not provided, return a zero-length numeric
   * vector. */
  if (length(x) == 0 || length(y) == 0)
    return allocVector(REALSXP, 0);

  /* filter for NULL and empty strings to make it easy to interface with R. */
  if (length(sx) == 0 || sx == R_NilValue)
    return utest(x, y, data, test, B, alpha, learning);
  else
    return ctest(x, y, sx, data, test, B, alpha, learning);

}/*INDEP_TEST*/
