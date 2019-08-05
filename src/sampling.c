#include "include/rcore.h"
#include "include/matrix.h"

void rbn_discrete_root(SEXP result, int cur, SEXP cpt, int *num, int ordinal,
    SEXP fixed);
void rbn_discrete_cond(SEXP result, SEXP nodes, int cur, SEXP parents, SEXP cpt,
    int *num, int ordinal, SEXP fixed);
void rbn_gaussian(SEXP result, int cur, SEXP parents, SEXP coefs, SEXP sigma,
    int *num, SEXP fixed);

/* sampling without replacement, internal copy of the SampleNoReplace function
 * in src/main/random.c.
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 1997--2010  The R Development Core Team
 * Copyright (C) 2003--2008  The R Foundation
 * licensed under "GPLv2 or later" licence. */
void SampleNoReplace(int k, int n, int *y, int *x) {

int i = 0, j = 0;

  for (i = 0; i < n; i++)
    x[i] = i;

  for (i = 0; i < k; i++) {

    j = n * unif_rand();
    y[i] = x[j] + 1;
    x[j] = x[--n];

  }/*FOR*/

}/*SAMPLENOREPLACE*/

/* sampling with replacement and equal probabilties, internal copy of the
 * SampleReplace function in src/main/random.c.
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 1997--2010  The R Development Core Team
 * Copyright (C) 2003--2008  The R Foundation
 * licensed under "GPLv2 or later" licence. */
void SampleReplace(int k, int n, int *y, int *x) {

int i;

  for (i = 0; i < k; i++)
    y[i] = x[(int)(n * unif_rand())];

}/*SAMPLEREPLACE*/

/* sampling with replacement and unequal probabilties, internal copy of the
 * ProbSampleReplace function in src/main/random.c.
 * Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 * Copyright (C) 1997--2010  The R Development Core Team
 * Copyright (C) 2003--2008  The R Foundation
 * licensed under "GPLv2 or later" licence. */
void ProbSampleReplace(int nprobs, double *probs, int *values, int ns,
    int *samples) {

double rU = 0;
int i = 0, j = 0;

  /* record element identities. */
  for (i = 0; i < nprobs; i++)
    values[i] = i + 1;

  /* sort the probabilities in descending order. */
  revsort(probs, values, nprobs);

  /* compute cumulative probabilities. */
  for (i = 1 ; i < nprobs; i++)
    probs[i] += probs[i - 1];

  /* generate the sample. */
  for (i = 0; i < ns; i++) {

    rU = unif_rand();

    for (j = 0; j < nprobs - 1; j++)
      if (rU <= probs[j])
        break;

    samples[i] = values[j];

  }/*FOR*/

}/*PROBSAMPLEREPLACE*/

void CondProbSampleReplace(int nprobs, int nconf, double *probs, int *conf,
    int *values, int ns, int *samples, bool *warn) {

int i = 0, j = 0;
double rU = 0;
bool *prepd = NULL;

  prepd = Calloc1D(nconf, sizeof(bool));

  /* generate the sample. */
  for (i = 0; i < ns; i++) {

    /* check whether the parents' configuration is missing. */
    if (conf[i] == NA_INTEGER) {

      samples[i] = NA_INTEGER;
      *warn = TRUE;
      continue;

    }/*THEN*/

    /* prepare the cumulative probabilities only for the conditional
       distributions we are effectively sampling from. */
    if (!prepd[conf[i]]) {

      /* sort the probabilities in descending order. */
      for (j = 0; j < nprobs; j++)
        values[CMC(j, conf[i], nprobs)] = j + 1;
      revsort(probs + conf[i] * nprobs, values + conf[i] * nprobs, nprobs);
      /* compute cumulative probabilities. */
      for (j = 1 ; j < nprobs; j++)
        probs[CMC(j, conf[i], nprobs)] += probs[CMC(j - 1, conf[i], nprobs)];
      /* flag the conditional distribution to avoid preparing it again. */
      prepd[conf[i]] = TRUE;

    }/*FOR*/

    /* check whether the conditional distribution is missing. */
    if (ISNAN(probs[CMC(0, conf[i], nprobs)])) {

      samples[i] = NA_INTEGER;
      *warn = TRUE;
      continue;

    }/*THEN*/

    rU = unif_rand();

    for (j = 0; j < nprobs; j++)
      if (rU <= probs[CMC(j, conf[i], nprobs)])
        break;

    samples[i] = values[CMC(j, conf[i], nprobs)];

  }/*FOR*/

  Free1D(prepd);

}/*CONDPROBSAMPLEREPLACE*/

