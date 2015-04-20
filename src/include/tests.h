
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

/* from mutual.information.c */
#define MI_PART(cell, xmarg, ymarg, zmarg) \
  ((cell) == 0 ? 0 : \
    ((double)(cell)) * log(((double)(cell)) * ((double)(zmarg)) / \
    (((double)(xmarg)) * ((double)(ymarg)))))

double c_mi(int *xx, int llx, int *yy, int lly, int num, double *df, int adj);
double c_cmi(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num,
    double *df, int adj);
double c_mig(double *xx, double *yy, int *num);
double c_micg(double *yy, double ym, double ysd, int *xx, int llx, int num);
double c_cmicg(double *yy, double **xx, int nx, int **zz, int nz, int *z0,
    int nz0, int *nlvls, int num);
double c_cmicg_unroll(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    double **gp, int ngp, double *df, int num);

/* from df.adjust.c */
double df_adjust(int *ni, int llx, int *nj, int lly);
double cdf_adjust(int **ni, int llx, int **nj, int lly, int llz);

/* from pearson.x2.c */
double c_x2(int *xx, int llx, int *yy, int lly, int num, double *df, int adj);
double c_cx2(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num,
    double *df, int adj);

/* from jonckheere.c */
double c_jt_mean(int num, int *ni, int llx);
double c_jt_var(int num, int *ni, int llx, int *nj, int lly);
double c_jt(int *xx, int llx, int *yy, int lly, int num);
double c_cjt(int *xx, int llx, int *yy, int lly, int *zz, int llz, int num);

/* from shrinkage.c */
void mi_lambda(double *n, double *lambda, double target, int num, int llx,
    int lly, int llz);
double c_shmi(int *xx, int llx, int *yy, int lly, int num);
double c_shcmi(int *xx, int llx, int *yy, int lly, int *zz, int llz,
    int num, double *df);
void c_cov_lambda(double **column, double *mean, int ncols, int n, double *var);
double c_fast_shcor(double *xx, double *yy, int *n);

/* from {discrete.gaussian}.monte.carlo.c */
#define MUTUAL_INFORMATION             1
#define PEARSON_X2                     2
#define SP_MUTUAL_INFORMATION          3
#define SP_PEARSON_X2                  4
#define JT                             5
#define GAUSSIAN_MUTUAL_INFORMATION    20
#define LINEAR_CORRELATION             21
#define FISHER_Z                       22
#define DISCRETE_PERMUTATION_TEST(t) t < GAUSSIAN_MUTUAL_INFORMATION

#define sequential_counter_check(counter) \
  (counter)++; \
  if ((counter) >= enough) { \
    (counter) = B; \
    break; \
  }

int remap_permutation_test(const char *t);
void c_mcarlo(int *xx, int nr, int *yy, int nc, int num, int B,
    double *observed, double *pvalue, double alpha, int test, double *df);
void c_cmcarlo(int *xx, int nr, int *yy, int nc, int *zz, int nl, int num,
    int B, double *observed, double *pvalue, double alpha, int test, double *df);
void c_gauss_mcarlo(double *xx, double *yy, int num, int B, double *res,
    double alpha, int test, double *observed);
void c_gauss_cmcarlo(double **column, int ncols, int num, int B,
    double *observed, double *pvalue, double alpha, int test);

