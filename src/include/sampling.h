
/* from sampling.c */
void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)
void SampleReplace(int k, int n, int *y, int *x);
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans);
void CondProbSampleReplace(int r, int c, double *p, int *conf, int *perm,
    int nans, int *ans, int *warn);

/* used to be in R core... */
void c_rcont2(int nrow, int ncol, int *nrowt, int *ncolt, int ntotal,
    double *fact, int *jwork, int **matrix);

/* from rbn.c */
void c_rbn_master(SEXP fitted, SEXP result, SEXP n, SEXP fix, int debuglevel);

/* from likelihood.weighting.c */
void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    int debuglevel);
