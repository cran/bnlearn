
/* from sampling.c */
void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)
void SampleReplace(int k, int n, int *y, int *x);
void ProbSampleReplace(int n, double *probs, int *values, int ns, int *samples);
void CondProbSampleReplace(int nprobs, int nconf, double *probs, int *conf,
    int *values, int ns, int *samples, bool *warn);

/* used to be in R core... */
void c_rcont2(int nrow, int ncol, int *nrowt, int *ncolt, int ntotal,
    double *fact, int *jwork, int **matrix);

/* from rbn.c */
void c_rbn_master(SEXP fitted, SEXP result, SEXP n, SEXP fix, bool debugging);

/* from likelihood.weighting.c */
void c_lw_weights(SEXP fitted, SEXP data, int n, double *w, SEXP keep,
    bool debugging);
