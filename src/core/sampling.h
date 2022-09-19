#ifndef SAMPLING_HEADER
#define SAMPLING_HEADER

void SampleNoReplace(int k, int n, int *y, int *x);
#define RandomPermutation(n, y, x) SampleNoReplace(n, n, y, x)
void SampleReplace(int k, int n, int *y, int *x);
void ProbSampleReplace(int n, double *probs, int *values, int ns, int *samples);
void CondProbSampleReplace(int nprobs, int nconf, double *probs, int *conf,
    int *values, int ns, int *samples, bool *warn);

#endif
