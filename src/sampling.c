#include "common.h"

void rbn_discrete_root(SEXP result, int cur, SEXP cpt, int *num);
void rbn_discrete_cond(SEXP result, int cur, SEXP cpt, int *configurations,
    int *num, int *warn);

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

    /* handle configurations which are not observed in the original data. */
    if (ISNAN(p[CMC(0, conf[i], r)]) || (conf[i] == NA_INTEGER)) {

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

/* generate random observations from a discrete bayesian network. */
SEXP rbn_discrete(SEXP fitted, SEXP n, SEXP debug) {

int *num = INTEGER(n), *poset = NULL, *debuglevel = LOGICAL(debug);
int *configurations = NULL;
int i = 0, k = 0, cur = 0, nnodes = LENGTH(fitted), nparents = 0, warn = 0;
SEXP result, nodes, roots, node_depth, cpt, parents, parent_vars, false;

  /* set up a logical variable to be used in the following calls. */
  PROTECT(false = allocVector(LGLSXP, 1));
  LOGICAL(false)[0] = FALSE;

  /* allocate and initialize the return value. */
  PROTECT(result = allocVector(VECSXP, nnodes));
  nodes = getAttrib(fitted, R_NamesSymbol);
  setAttrib(result, R_NamesSymbol, nodes);

  /* order the nodes according to their depth in the graph. */
  PROTECT(roots = root_nodes(fitted, false));
  PROTECT(node_depth = schedule(fitted, roots, false, false));
  poset = alloc1dcont(nnodes);
  for (i = 0; i < nnodes; i++)
    poset[i] = i;
  R_qsort_int_I(INTEGER(node_depth), poset, 1, nnodes);

  /* unprotect roots and node_depth, they are not needed any more. */
  UNPROTECT(2);

  if (*debuglevel > 0) {

    Rprintf("* partial node ordering is:");

    for (i = 0; i < nnodes; i++)
      Rprintf(" %s", NODE(poset[i]));

    Rprintf(".\n");

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  /* allocate the working space for the configurations' indexes. */
  configurations = alloc1dcont(*num);

  for (i = 0; i < nnodes; i++) {

    /* get the index of the node we have to generate random observations from,
     * its conditional probability table and the number of its parents. */
    cur = poset[i];
    cpt = getListElement(VECTOR_ELT(fitted, cur), "prob");
    parents = getListElement(VECTOR_ELT(fitted, cur), "parents");
    nparents = LENGTH(parents);

    /* generate the random observations for the current node. */
    if (nparents == 0) {

      if (*debuglevel > 0)
        Rprintf("* simulating node %s, which doesn't have any parent.\n", NODE(cur));

      rbn_discrete_root(result, cur, cpt, num);

    }/*THEN*/
    else {

      if (*debuglevel > 0) {

        Rprintf("* simulating node %s with parents ", NODE(cur));
        for (k = 0; k < nparents - 1; k++)
          Rprintf("%s, ", CHAR(STRING_ELT(parents, k)));
        Rprintf("%s.\n", CHAR(STRING_ELT(parents, nparents - 1)));

      }/*THEN*/

      PROTECT(parent_vars = dataframe_column(result, parents, false));
      cfg(parent_vars, configurations, NULL);

      rbn_discrete_cond(result, cur, cpt, configurations, num, &warn);

      if (warn == TRUE) {

        warning("some configurations of the parents of %s are not present in the original data. NAs will be generated.",
          NODE(cur));

      }/*THEN*/

      UNPROTECT(1);

    }/*ELSE*/

  }/*FOR*/

  PutRNGstate();

  /* add the labels to the return value. */
  minimal_data_frame(result);

  UNPROTECT(2);
  return result;

}/*RBN_DISCRETE*/

/* this is a root node, so this a classic, unconditional sampling. */
void rbn_discrete_root(SEXP result, int cur, SEXP cpt, int *num) {

int np = LENGTH(cpt), *workplace = NULL;
double *p = NULL;
SEXP generated, class, lvls;

  workplace = alloc1dcont(np);

  /* allocate the memory for the generated observations. */
  PROTECT(generated = allocVector(INTSXP, *num));
  /* duplicate the probability table to save the original copy from tampering. */
  p = alloc1dreal(np);
  memcpy(p, REAL(cpt), np * sizeof(double));

  /* perform the random sampling. */
  ProbSampleReplace(np, p, workplace, *num, INTEGER(generated));
  /* set up all the attributes for the newly generated observations. */
  lvls = VECTOR_ELT(getAttrib(cpt, R_DimNamesSymbol), 0);
  setAttrib(generated, R_LevelsSymbol, lvls);
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(generated, R_ClassSymbol, class);

  /* add the resulting vector to the return value of rbn_discrete(). */
  SET_VECTOR_ELT(result, cur, generated);

  UNPROTECT(2);

}/*RBN_DISCRETE_ROOT*/

/* this is not a root node, so this is a conditional sampling. */
void rbn_discrete_cond(SEXP result, int cur, SEXP cpt, int *configurations,
    int *num, int *warn) {

int np = LENGTH(cpt), nlevels = 0, *workplace = NULL;
double *p = NULL;
SEXP generated, class, lvls;

  /* get the number of levels of the curent variable .*/
  lvls = VECTOR_ELT(getAttrib(cpt, R_DimNamesSymbol), 0);
  nlevels = LENGTH(lvls);

  workplace = alloc1dcont(np);

  /* allocate the memory for the generated observations. */
  PROTECT(generated = allocVector(INTSXP, *num));
  /* duplicate the probability table to save the original copy from tampering. */
  p = alloc1dreal(np);
  memcpy(p, REAL(cpt), np * sizeof(double));
  /* perform the random sampling. */
  CondProbSampleReplace(nlevels, LENGTH(cpt)/nlevels, p, configurations,
    workplace, *num, INTEGER(generated), warn);
  /* set up all the attributes for the newly generated observations. */
  setAttrib(generated, R_LevelsSymbol, lvls);
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("factor"));
  setAttrib(generated, R_ClassSymbol, class);

  /* add the resulting vector to the return value of rbn_discrete(). */
  SET_VECTOR_ELT(result, cur, generated);

  UNPROTECT(2);

}/*RBN_DISCRETE_COND*/

