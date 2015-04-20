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
void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans) {

double rU = 0;
int i = 0, j = 0;
int nm1 = n - 1;

  /* record element identities. */
  for (i = 0; i < n; i++)
    perm[i] = i + 1;

  /* sort the probabilities into descending order. */
  revsort(p, perm, n);

  /* compute cumulative probabilities. */
  for (i = 1 ; i < n; i++)
    p[i] += p[i - 1];

  /* compute the sample. */
  for (i = 0; i < nans; i++) {

    rU = unif_rand();

    for (j = 0; j < nm1; j++)
      if (rU <= p[j])
        break;

    ans[i] = perm[j];

  }/*FOR*/

}/*PROBSAMPLEREPLACE*/

void CondProbSampleReplace(int r, int c, double *p, int *conf, int *perm,
    int nans, int *ans, int *warn) {

int i = 0, j = 0;
double rU = 0;

  /* record element identities. */
  for (i = 0; i < r; i ++)
    for (j = 0; j < c; j++)
      perm[CMC(i, j, r)] = i + 1;

  /* sort the probabilities into descending order. */
  for (j = 0; j < c; j++)
    revsort(p + j * r, perm + j * r, r);

  /* compute cumulative probabilities. */
  for (j = 0; j < c; j++)
    for (i = 1 ; i < r; i++)
      p[CMC(i, j, r)] += p[CMC(i - 1, j, r)];

  /* compute the sample. */
  for (i = 0; i < nans; i++) {

    /* check whether the parents' configuration is missing. */
    if (conf[i] == NA_INTEGER) {

      ans[i] = NA_INTEGER;
      *warn = TRUE;
      continue;

    }/*THEN*/

    /* check whether the conditional distribution is missing. */
    if (ISNAN(p[CMC(0, conf[i], r)])) {

      ans[i] = NA_INTEGER;
      *warn = TRUE;
      continue;

    }/*THEN*/

    rU = unif_rand();

    for (j = 0; j < r; j++)
      if (rU <= p[CMC(j, conf[i], r)])
        break;

    ans[i] = perm[CMC(j, conf[i], r)];

  }/*FOR*/

}/*CONDPROBSAMPLEREPLACE*/

