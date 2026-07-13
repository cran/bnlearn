#ifndef TEST_PATTERNS_HEADER
#define TEST_PATTERNS_HEADER

#include "../core/data.table.h"
#include "tests.h"

/* unconditional tests. */
double ut_discrete(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df, test_e test);
double ut_gaustests_complete(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test);
double ut_gaustests_with_missing(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df, test_e test);
double ut_micg_complete(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df);
double ut_micg_with_missing(SEXP xx, SEXP yy, int nobs, int ntests,
    double *pvalue, double *df);
double ut_dperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    double *df, test_e type, int B, double a);
double ut_gperm(SEXP xx, SEXP yy, int nobs, int ntests, double *pvalue,
    test_e type, int B, double a, bool complete);
double ut_custom(SEXP x, SEXP y, SEXP data, SEXP custom_fn, SEXP custom_args,
    double *pvalue);

/* conditional tests. */
double ct_discrete(tabular dtx, tabular dty, tabular dtz, double *pvalue, double *df,
    test_e test);
double ct_gaustests_complete(tabular dtx, tabular dt, double *pvalue, double *df,
    test_e test);
double ct_gaustests_with_missing(tabular dtx, tabular dt, double *pvalue,
    double *df, test_e test);
double c_micg_chisq(tabular dtx, tabular dty, tabular dtz, int *zptr, int llz,
    double *df, bool copy);
double ct_micg_complete(tabular dtx, tabular dty, tabular dtz, double *pvalue,
    double *df);
double ct_micg_with_missing(tabular dtx, tabular dty, tabular dt, double *pvalue,
    double *df);
double ct_dperm(tabular dtx, tabular dty, tabular dtz, double *pvalue, double *df,
    test_e type, int B, double a);
double ct_gperm(tabular dtx, tabular dt, double *pvalue, double *df, test_e type,
    int B, double a, bool complete);
double ct_custom(SEXP x, SEXP y, SEXP z, SEXP data, SEXP custom_fn,
    SEXP custom_args, double *pvalue);

/* round-robin tests. */
void rrd_discrete(tabular dtx, tabular dtz, test_e test, double *pvalue,
    double alpha, bool debugging);
void rrd_gaustests_complete(tabular dt, test_e test, double *pvalue, double
    alpha, bool debugging);
void rrd_gaustests_with_missing(tabular dt, test_e test, double *pvalue,
    double alpha, bool debugging);
double rrd_micg_chisq(tabular dtx, tabular dty, tabular sub, int *zptr, int llz,
    double *df, bool copy);
void rrd_micg_complete(tabular dtx, tabular dtz, test_e test, double *pvalue,
    double alpha, bool debugging);
void rrd_micg_with_missing(tabular dtx, tabular dtz, test_e test, double *pvalue,
    double alpha, bool debugging);
void rrd_dperm(tabular dtx, tabular dtz, test_e test, double *pvalue, double alpha,
    int nperms, double threshold, bool debugging);
void rrd_gperm(tabular dt, test_e test, double *pvalue, double alpha, int nperms,
    double threshold, bool complete, bool debugging);
void rrd_custom(SEXP x, SEXP z, SEXP fixed, SEXP data, SEXP custom_fn,
    SEXP custom_args, double *pvalue, double alpha, bool debugging);

/* all-subsets tests. */
SEXP ast_discrete(tabular dtx, tabular dty, tabular dtz, int nf, int minsize,
    int maxsize, test_e test, double a, bool debugging);
SEXP ast_gaustests_complete(tabular dt, int nf, int minsize, int maxsize,
    double a, bool debugging, test_e test);
SEXP ast_gaustests_with_missing(tabular dt, int nf, int minsize, int maxsize,
    double a, bool debugging, test_e test);
SEXP ast_micg_complete(tabular dtx, tabular dty, tabular dtz, int nf, int minsize,
   int maxsize, double a, bool debugging);
SEXP ast_micg_with_missing(tabular dtx, tabular dty, tabular dtz, int nf,
    int minsize, int maxsize, double a, bool debugging);
SEXP ast_dperm(tabular dtx, tabular dty, tabular dtz, int nf, int minsize,
    int maxsize, double a, test_e type, int nperms, double threshold,
    bool debugging);
SEXP ast_gperm(tabular dt, int nf, int minsize, int maxsize, double a,
    test_e type, int nperms, double threshold, bool complete, bool debugging);
SEXP ast_custom(SEXP x, SEXP y, SEXP sx, SEXP fixed, SEXP data, int minsize,
    int maxsize, double a, SEXP custom_fn, SEXP custom_args,
    bool debugging);

SEXP ast_prepare_retval(double pvalue, double min_pvalue, double max_pvalue,
    double alpha, const char **nodes, int nnodes);
void update_pvalue_range(double pvalue, double *min, double *max);

#endif
