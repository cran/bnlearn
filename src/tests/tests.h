#ifndef TEST_STATISTICS_HEADER
#define TEST_STATISTICS_HEADER

#include "../core/covariance.matrix.h"
#include "../core/contingency.tables.h"

/* enum for tests, to be matched from the label string passed down from R. */
typedef enum {
  ENOTEST   =  0, /* error code, no such test. */
  MI        =  1, /* mutual information test, without df adjustments. */
  MI_ADF    =  2, /* mutual information test, with df adjustments. */
  X2        =  3, /* Pearson's X^2 test, without df adjustments. */
  X2_ADF    =  4, /* Pearson's X^2 test, with df adjustments. */
  JT        = 10, /* Jonckheere-Terpstra test. */
  COR       = 20, /* linear correlation. */
  ZF        = 21, /* Fisher's Z test. */
  MI_G      = 22, /* Gaussian mutual information test. */
  MI_CG     = 30, /* conditional linear Gaussian mutual information test. */
  MI_SH     = 40, /* shrinkage mutual information test. */
  MI_G_SH   = 41, /* shrinkage Gaussian mutual information test. */
  MC_MI     = 50, /* mutual information test (permutation). */
  MC_X2     = 51, /* Pearson's X^2 test (permutation). */
  SP_MI     = 52, /* mutual information test (semiparametric). */
  SP_X2     = 53, /* Pearson's X^2 test (semiparametric). */
  MC_JT     = 54, /* Jonckheere-Terpstra test (permutation). */
  SMC_MI    = 60, /* mutual information test (sequential permutation). */
  SMC_X2    = 61, /* Pearson's X^2 test (sequential permutation). */
  SMC_JT    = 62, /* Jonckheere-Terpstra test (sequential permutation). */
  MC_COR    = 70, /* linear correlation (permutation). */
  MC_MI_G   = 71, /* Gaussian mutual information test (permutation). */
  MC_ZF     = 72, /* Fisher's Z test (permutation). */
  SMC_COR   = 80, /* linear correlation (sequential permutation). */
  SMC_MI_G  = 81, /* Gaussian mutual information test (sequential permutation). */
  SMC_ZF    = 82, /* Fisher's Z test (sequential permutation). */
  CUSTOM_T  = 99, /* user-provided test function. */
} test_e;

/* from enums.c */
test_e test_to_enum(const char *label);
#define IS_DISCRETE_ASYMPTOTIC_TEST(t) (((t >= MI) && (t < COR)) || t == MI_SH)
#define IS_DISCRETE_PERMUTATION_TEST(t) ((t >= MC_MI) && (t < MC_COR))
#define IS_CONTINUOUS_PERMUTATION_TEST(t) ((t >= MC_COR) && (t < CUSTOM_T))
#define IS_SMC(t) (((t >= SMC_MI) && (t < MC_COR)) || ((t >= SMC_COR) && (t < CUSTOM_T)))
#define IS_TWO_SIDED(t) \
  ((t == JT) || (t == COR) || (t == ZF) || (t == MC_JT) || (t == SMC_JT) || \
   (t == MC_COR) || (t == MC_ZF) || (t == SMC_COR) || (t == SMC_ZF))

#define SEQUENTIAL_COUNTER_CHECK(counter) \
  { \
    (counter)++; \
    if ((counter) >= enough) { \
	  (counter) = B; \
	  break; \
    } \
  }

/* from htest.c */
SEXP c_create_htest(double stat, SEXP test, double pvalue, double df,
    SEXP extra_args);

/* from indep.test.c */
SEXP indep_test(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP learning, SEXP complete);

/* conditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP alpha, SEXP extra_args,
    SEXP learning, SEXP complete);
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP alpha,
    SEXP extra_args, SEXP learning, SEXP complete);

/* from discrete.tests.c */
double c_chisqtest(int *xx, int llx, int *yy, int lly, int num, double *df,
    test_e test, bool scale);
double mi_kernel(counts2d table);
double x2_kernel(counts2d table);
double mi_kernel_collapsed(counts2d table, int k);
double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test, bool scale);
double cmi_kernel(counts3d table);
double cx2_kernel(counts3d table);

/* from gaussian.tests.c */
double cor_t_trans(double cor, double df);
double cor_zf_trans(double cor, double df);
double cor_mi_trans(double cor);

/* from cg.mutual.information.c */
double c_micg(double *yy, double ym, double ysd, int *xx, int llx, int num,
    double *df);
double c_cmicg(double *yy, double **xx, int nx, int **zz, int nz, int *z0,
    int nz0, int *nlvls, int num, double *df);
double c_cmicg_unroll(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    double **gp, int ngp, double *df, int num);
double c_micg_with_missing(double *yy, int *xx, int llx, int num, double *df,
    int *ncomplete);

/* from df.adjust.c */
double discrete_df(test_e test, int *ni, int llx, int *nj, int lly);
double discrete_cdf(test_e test, int **ni, int llx, int **nj, int lly, int llz);
double gaussian_cdf(test_e test, int num, int nz);
#define gaussian_df(test, num) gaussian_cdf(test, num, 0)

/* from jonckheere.c */
double c_jt(int *xx, int llx, int *yy, int lly, int num);
double c_cjt(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num);
double jt_centered_kernel(counts2d table);
double cjt_centered_kernel(counts3d table);
double jt_var_kernel(counts2d table);
double cjt_var_kernel(counts3d table);

/* from shrinkage.c */
void mi_lambda(double *n, double *lambda, double target, int num, int llx,
    int lly, int llz);
double cor_lambda(double *xx, double *yy, int nobs, int ncomplete,
   double xm, double ym, double xsd, double ysd, double cor);
double c_shmi(int *xx, int llx, int *yy, int lly, int num, int scale);
double c_shcmi(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, int scale);
double covmat_lambda(double **column, double *mean, covariance cov, int n,
    bool *missing, int nc);
void covmat_shrink(covariance cov, double lambda);

/* from {discrete,gaussian}.monte.carlo.c */
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, test_e test, double *df);
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, test_e test, double *df);
void c_gauss_mcarlo(double *xx, double *yy, int num, int B, double *res,
    double alpha, test_e test, double *observed);
void c_gauss_cmcarlo(double **column, int ncol, int num, int v1, int v2, int B,
    double *observed, double *pvalue, double alpha, test_e test);

/* from custom.test.c */
double custom_test_function(SEXP x, SEXP y, SEXP z, SEXP data, SEXP custom_fn,
    SEXP custom_args, double *pvalue);

#endif
