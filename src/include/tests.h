
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
} test_e;

/* from enums.c */
test_e test_label(const char *label);
#define IS_DISCRETE_ASYMPTOTIC_TEST(t) (((t >= MI) && (t < COR)) || t == MI_SH)
#define IS_DISCRETE_PERMUTATION_TEST(t) ((t >= MC_MI) && (t < MC_COR))
#define IS_CONTINUOUS_PERMUTATION_TEST(t) (t >= MC_COR)
#define IS_SMC(t) (((t >= SMC_MI) && (t < MC_COR)) || (t >= SMC_COR))
#define IS_ADF(t) ((t == MI_ADF) || (t == X2_ADF))
#define IS_TWO_SIDED(t) \
  ((t == JT) || (t == COR) || (t == ZF) || (t == MC_JT) || (t == SMC_JT) || \
   (t == MC_COR) || (t == MC_ZF) || (t == SMC_COR) || (t == SMC_ZF))

#define SEQUENTIAL_COUNTER_CHECK(counter) \
  (counter)++; \
  if ((counter) >= enough) { \
    (counter) = B; \
    break; \
  }

/* from htest.c */
SEXP c_create_htest(double stat, SEXP test, double pvalue, double df, SEXP B);

/* from indep.test.c */
SEXP indep_test(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B,
    SEXP alpha, SEXP learning);

/* conditional independence tests. */
SEXP utest(SEXP x, SEXP y, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning);
SEXP ctest(SEXP x, SEXP y, SEXP sx, SEXP data, SEXP test, SEXP B, SEXP alpha,
    SEXP learning);

/* from discrete.tests.c */
double c_chisqtest(int *xx, int llx, int *yy, int lly, int num, double *df,
    test_e test);
double mi_kernel(int **n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length);
double x2_kernel(int **n, int *nrowt, int *ncolt, int nrows, int ncols,
    int length);
double c_cchisqtest(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df, test_e test);
double cmi_kernel(int ***n, int **nrowt, int **ncolt, int *ncond, int nr,
    int nc, int nl);
double cx2_kernel(int ***n, int **nrowt, int **ncolt, int *ncond,
    int nr, int nc, int nl);

/* from gaussian.tests.c */
double cor_t_trans(double cor, double df);
double cor_zf_trans(double cor, double df);
double cor_mi_trans(double cor);

/* from cg.mutual.information.c */
double c_micg(double *yy, double ym, double ysd, int *xx, int llx, int num);
double c_cmicg(double *yy, double **xx, int nx, int **zz, int nz, int *z0,
    int nz0, int *nlvls, int num);
double c_cmicg_unroll(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    double **gp, int ngp, double *df, int num);

/* from df.adjust.c */
double df_adjust(int *ni, int llx, int *nj, int lly);
double cdf_adjust(int **ni, int llx, int **nj, int lly, int llz);

/* from jonckheere.c */
double c_jt_mean(int num, int *ni, int llx);
double c_jt_var(int num, int *ni, int llx, int *nj, int lly);
double c_jt(int *xx, int llx, int *yy, int lly, int num);
double c_cjt(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num);

/* from shrinkage.c */
void mi_lambda(double *n, double *lambda, double target, int num, int llx,
    int lly, int llz);
double cor_lambda(double *xx, double *yy, int nobs, double xm, double ym,
   double xsd, double ysd, double cor);
double c_shmi(int *xx, int llx, int *yy, int lly, int num);
double c_shcmi(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df);
double covmat_lambda(double **column, double *mean, double *var, int n,
    int ncols);
void covmat_shrink(double *var, int ncols, double lambda);

/* from {discrete,gaussian}.monte.carlo.c */
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, test_e test, double *df);
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, test_e test, double *df);
void c_gauss_mcarlo(double *xx, double *yy, int num, int B, double *res,
    double alpha, test_e test, double *observed);
void c_gauss_cmcarlo(double **column, int ncols, int num, int B,
    double *observed, double *pvalue, double alpha, test_e test);

/* from contingency.tables.c */
void fill_1d_table(int *xx, int **n, int llx, int num);
void fill_2d_table(int *xx, int *yy, int ***n, int **ni, int **nj, int llx,
    int lly, int num);
void fill_3d_table(int *xx, int *yy, int *zz, int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz, int num);

