#include "common.h"

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

#define CATEGORICAL      1
#define ORDINAL          2
#define GAUSSIAN         3

SEXP rbn_expand_fix(SEXP fix, SEXP nodes, int *nnodes) {

int i = 0, *f = NULL;
SEXP result, fixed_nodes, try;

  /* get the names of the nodes that are fixed. */
  fixed_nodes = getAttrib(fix, R_NamesSymbol);
  /* match the names with the nodes in the network. */
  PROTECT(try = match(nodes, fixed_nodes, 0));
  f = INTEGER(try);
  /* allocate the return value. */
  PROTECT(result = allocVector(VECSXP, *nnodes));
  setAttrib(result, R_NamesSymbol, nodes);

  for (i = 0; i < LENGTH(fixed_nodes); i++)
    SET_VECTOR_ELT(result, f[i] - 1, VECTOR_ELT(fix, i));

  UNPROTECT(2);
  return result;

}/*RBN_EXPAND_FIX*/

/* generate random observations from a bayesian network. */
SEXP rbn_master(SEXP fitted, SEXP n, SEXP fix, SEXP debug) {

int *num = INTEGER(n), *poset = NULL, *debuglevel = LOGICAL(debug), type = 0;
int has_fixed = (TYPEOF(fix) != LGLSXP);
int i = 0, k = 0, cur = 0, nnodes = LENGTH(fitted), nparents = 0;
const char *cur_class = NULL;
SEXP result, nodes, roots, node_depth, cpt, coefs, sd, parents, parent_vars, false;
SEXP cur_node, cur_fixed;

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

  if (has_fixed)
    PROTECT(fix = rbn_expand_fix(fix, nodes, &nnodes));

  if (*debuglevel > 0) {

    Rprintf("* partial node ordering is:");

    for (i = 0; i < nnodes; i++)
      Rprintf(" %s", NODE(poset[i]));

    Rprintf(".\n");

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  for (i = 0; i < nnodes; i++) {

    /* get the index of the node we have to generate random observations from,
     * its conditional probability table/regression parameters and the number
     * of its parents. */
    cur = poset[i];
    cur_node = VECTOR_ELT(fitted, cur);
    cur_class = CHAR(STRING_ELT(getAttrib(cur_node, R_ClassSymbol), 0));
    parents = getListElement(cur_node, "parents");
    nparents = LENGTH(parents);

    /* check whether the value of the node is fixed, and if so retrieve it from 
     * the list. */
    if (has_fixed)
      cur_fixed = VECTOR_ELT(fix, cur);
    else
      cur_fixed = R_NilValue;

    /* find out whether the node corresponds to an ordered factor or not. */
    if (strcmp(cur_class, "bn.fit.onode") == 0) {

      cpt = getListElement(cur_node, "prob");
      type = ORDINAL;

    }/*THEN*/
    else if (strcmp(cur_class, "bn.fit.dnode") == 0) {

      cpt = getListElement(cur_node, "prob");
      type = CATEGORICAL;

    }/*THEN*/
    else if (strcmp(cur_class, "bn.fit.gnode") == 0) {

      coefs = getListElement(cur_node, "coefficients");
      sd = getListElement(cur_node, "sd");
      type = GAUSSIAN;

    }/*THEN*/

    /* generate the random observations for the current node. */
    if (nparents == 0) {

      if (*debuglevel > 0) {

        if (cur_fixed != R_NilValue)
          Rprintf("* node %s is fixed.\n", NODE(cur));
        else
          Rprintf("* simulating node %s, which doesn't have any parent.\n",
            NODE(cur));

      }/*THEN*/

      switch(type) {

        case CATEGORICAL:
          rbn_discrete_root(result, cur, cpt, num, FALSE, cur_fixed);
          break;

        case ORDINAL:
          rbn_discrete_root(result, cur, cpt, num, TRUE, cur_fixed);
          break;

        case GAUSSIAN:
          rbn_gaussian(result, cur, NULL, coefs, sd, num, cur_fixed);
          break;

      }/*SWITCH*/

    }/*THEN*/
    else {

      if (*debuglevel > 0) {

        if (cur_fixed != R_NilValue) {

          Rprintf("* node %s is fixed, ignoring parents.\n", NODE(cur));

        }/*THEN*/
        else {

          Rprintf("* simulating node %s with parents ", NODE(cur));
          for (k = 0; k < nparents - 1; k++)
            Rprintf("%s, ", CHAR(STRING_ELT(parents, k)));
          Rprintf("%s.\n", CHAR(STRING_ELT(parents, nparents - 1)));

        }/*ELSE*/

      }/*THEN*/

      PROTECT(parent_vars = dataframe_column(result, parents, false));

      switch(type) {

        case CATEGORICAL:
          rbn_discrete_cond(result, nodes, cur, parent_vars, cpt, num, FALSE, cur_fixed);
          break;

        case ORDINAL:
          rbn_discrete_cond(result, nodes, cur, parent_vars, cpt, num, TRUE, cur_fixed);
          break;

        case GAUSSIAN:
          rbn_gaussian(result, cur, parent_vars, coefs, sd, num, cur_fixed);
          break;

      }/*SWITCH*/

      UNPROTECT(1);

    }/*ELSE*/

  }/*FOR*/

  PutRNGstate();

  /* add the labels to the return value. */
  minimal_data_frame(result);

  UNPROTECT(2 + has_fixed);
  return result;

}/*RBN_MASTER*/

/* unconditional discrete sampling. */
void rbn_discrete_root(SEXP result, int cur, SEXP cpt, int *num, int ordinal,
    SEXP fixed) {

int np = LENGTH(cpt), *gen = NULL, *workplace = NULL;
double *p = NULL;
SEXP generated, class, lvls;

  /* get the levels of the curent variable .*/
  lvls = VECTOR_ELT(getAttrib(cpt, R_DimNamesSymbol), 0);
  /* allocate the memory for the generated observations. */
  PROTECT(generated = allocVector(INTSXP, *num));
  gen = INTEGER(generated);

  if (fixed != R_NilValue) {

    int constant = INTEGER(match(lvls, fixed, 0))[0];

    for (int i = 0; i < *num; i++)
      gen[i] = constant;

  }/*THEN*/
  else {

    workplace = alloc1dcont(np);

    /* duplicate the probability table to save the original copy from tampering. */
    p = alloc1dreal(np);
    memcpy(p, REAL(cpt), np * sizeof(double));

    /* perform the random sampling. */
    ProbSampleReplace(np, p, workplace, *num, gen);

  }/*ELSE*/

  /* set up all the attributes for the newly generated observations. */
  setAttrib(generated, R_LevelsSymbol, lvls);
  if (ordinal) {

    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("ordered"));
    SET_STRING_ELT(class, 1, mkChar("factor"));

  }/*THEN*/
  else {

    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("factor"));

  }/*ELSE*/
  setAttrib(generated, R_ClassSymbol, class);

  /* add the resulting vector to the return value of rbn_discrete(). */
  SET_VECTOR_ELT(result, cur, generated);

  UNPROTECT(2);

}/*RBN_DISCRETE_ROOT*/

/* conditional discrete sampling. */
void rbn_discrete_cond(SEXP result, SEXP nodes, int cur, SEXP parents, SEXP cpt,
    int *num, int ordinal, SEXP fixed) {

int np = LENGTH(cpt), nlevels = 0, warn = 0;
int *workplace = NULL, *configurations = NULL, *gen = NULL;
double *p = NULL;
SEXP generated, class, lvls;

  /* get the number of levels of the curent variable .*/
  lvls = VECTOR_ELT(getAttrib(cpt, R_DimNamesSymbol), 0);
  nlevels = LENGTH(lvls);
  /* allocate the memory for the generated observations. */
  PROTECT(generated = allocVector(INTSXP, *num));
  gen = INTEGER(generated);

  if (fixed != R_NilValue) {

    int constant = INTEGER(match(lvls, fixed, 0))[0];

    for (int i = 0; i < *num; i++)
      gen[i] = constant;

  }/*THEN*/
  else {

    workplace = alloc1dcont(np);

    /* allocate and initialize the parents' configurations. */
    configurations = alloc1dcont(*num);
    cfg(parents, configurations, NULL);

    /* duplicate the probability table to save the original copy from tampering. */
    p = alloc1dreal(np);
    memcpy(p, REAL(cpt), np * sizeof(double));
    /* perform the random sampling. */
    CondProbSampleReplace(nlevels, LENGTH(cpt)/nlevels, p, configurations,
      workplace, *num, gen, &warn);

  }/*ELSE*/

  /* set up all the attributes for the newly generated observations. */
  setAttrib(generated, R_LevelsSymbol, lvls);
  if (ordinal) {

    PROTECT(class = allocVector(STRSXP, 2));
    SET_STRING_ELT(class, 0, mkChar("ordered"));
    SET_STRING_ELT(class, 1, mkChar("factor"));

  }/*THEN*/
  else {

    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("factor"));

  }/*ELSE*/
  setAttrib(generated, R_ClassSymbol, class);

  /* add the resulting vector to the return value of rbn_discrete(). */
  SET_VECTOR_ELT(result, cur, generated);

  /* warn when returning missing values. */
  if (warn == TRUE)
    warning("some configurations of the parents of %s are not present in the original data. NAs will be generated.",
      NODE(cur));

  UNPROTECT(2);

}/*RBN_DISCRETE_COND*/

/* conditional and unconditional normal sampling. */
void rbn_gaussian(SEXP result, int cur, SEXP parents, SEXP coefs, SEXP sigma,
    int *num, SEXP fixed) {

int i = 0, j = 0, p = LENGTH(coefs);
double *beta = REAL(coefs), *sd = REAL(sigma), *gen = NULL, *Xj = NULL;
SEXP generated;

  /* allocate the memory for the generated observations. */
  PROTECT(generated = allocVector(REALSXP, *num));
  gen = REAL(generated);

  if (fixed != R_NilValue) {

    double *constant = REAL(fixed);

    for (i = 0; i < *num; i++)
      gen[i] = *constant;

  }/*THEN*/
  else {

    /* initialize with intercept and standard error. */
    for (i = 0; i < *num; i++)
      gen[i] = beta[0] + norm_rand() * (*sd);

    /* add the contributions of the other regressors (if any). */
    for (j = 1; j < p; j++) {

      Xj = REAL(VECTOR_ELT(parents, j - 1));

      for (i = 0; i < *num; i++) 
        gen[i] += Xj[i] * beta[j];

    }/*FOR*/

  }/*ELSE*/

  /* add the resulting vector to the return value of rbn_discrete(). */
  SET_VECTOR_ELT(result, cur, generated);

  UNPROTECT(1);

}/*RBN_GAUSSIAN*/

