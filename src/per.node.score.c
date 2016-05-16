#include "include/rcore.h"
#include "include/scores.h"

#define DEBUG_BEFORE() \
  if (debuglevel > 0) { \
    Rprintf("----------------------------------------------------------------\n"); \
    Rprintf("* processing node %s.\n", CHAR(STRING_ELT(cur, 0))); \
  }/*THEN*/

/* R frontend: compute the score component for each target node. */
SEXP per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, SEXP debug) {

SEXP result;

  /* allocate the return value. */
  PROTECT(result = allocVector(REALSXP, length(targets)));
  /* compute the score componenets. */
  c_per_node_score(network, data, score, targets, extra_args, isTRUE(debug),
    REAL(result));
  /* set labels on the computed score components. */
  setAttrib(result, R_NamesSymbol, targets);

  UNPROTECT(1);

  return result;

}/*PER_NODE_SCORE*/

/* C backend: compute the score component for each target node. */
void c_per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, int debuglevel, double *res) {

int i = 0, ntargets = length(targets);
char *s = (char *)CHAR(STRING_ELT(score, 0));
SEXP cur;

  /* allocate dummy variable for the current node's label. */
  PROTECT(cur = allocVector(STRSXP, 1));

  if (strcmp(s, "loglik") == 0) {

    /* discrete log-likelihood score. */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_dnode(cur, network, data, NULL, debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(s, "loglik-g") == 0) {

    /* Gaussian log-likelihood score. */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_gnode(cur, network, data, NULL, debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(s, "loglik-cg") == 0) {

    /* Conditional Linear Gaussian log-likelihood score. */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_cgnode(cur, network, data, NULL, debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(s, "aic") == 0) || (strcmp(s, "bic") == 0)) {

    /* AIC and BIC scores, discrete data. */
    double nparams = 0, *k = REAL(getListElement(extra_args, "k"));

    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_dnode(cur, network, data, &nparams, debuglevel);
      res[i] -= (*k) * nparams;

      if (debuglevel > 0)
        Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(s, "aic-g") == 0) || (strcmp(s, "bic-g") == 0)) {

    /* AIC and BIC scores, Gaussian data. */
    double nparams = 0, *k = REAL(getListElement(extra_args, "k"));

    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_gnode(cur, network, data, &nparams, debuglevel);
      res[i] -= (*k) * nparams;

      if (debuglevel > 0)
        Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(s, "aic-cg") == 0) || (strcmp(s, "bic-cg") == 0)) {

    /* AIC and BIC scores, Conditional Linear Gaussian data. */
    double nparams = 0, *k = REAL(getListElement(extra_args, "k"));

    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = loglik_cgnode(cur, network, data, &nparams, debuglevel);
      res[i] -= (*k) * nparams;

      if (debuglevel > 0)
        Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

    }/*FOR*/

  }/*THEN*/
  else if ((strcmp(s, "bde") == 0) || (strcmp(s, "bds") == 0)) {

    SEXP iss = getListElement(extra_args, "iss");
    SEXP prior = getListElement(extra_args, "prior");
    SEXP beta = getListElement(extra_args, "beta");

    /* Bayesian Dirichlet equivalent score (BDe) and sparse score (BDs). */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = dirichlet_node(cur, network, data, iss, prior, beta, R_NilValue,
                 (strcmp(s, "bds") == 0), debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(s, "k2") == 0) {

    /* K2 score. */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = dirichlet_node(cur, network, data, R_NilValue, R_NilValue, R_NilValue,
                 R_NilValue, FALSE, debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(s, "bge") == 0) {

    SEXP iss = getListElement(extra_args, "iss");
    SEXP phi = getListElement(extra_args, "phi");
    SEXP prior = getListElement(extra_args, "prior");
    SEXP beta = getListElement(extra_args, "beta");

    /* Bayesian Gaussian equivalent score (BGe). */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = wishart_node(cur, network, data, iss, phi, prior, beta, debuglevel);

    }/*FOR*/

  }/*THEN*/
  else if (strcmp(s, "mbde") == 0) {

    SEXP iss = getListElement(extra_args, "iss");
    SEXP exp = getListElement(extra_args, "exp");
    SEXP prior = getListElement(extra_args, "prior");
    SEXP beta = getListElement(extra_args, "beta");

    /* Mixture Bayesian Dirichlet equivalent score (mBDe). */
    for (i = 0; i < ntargets; i++) {

      SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
      DEBUG_BEFORE();
      res[i] = dirichlet_node(cur, network, data, iss, prior, beta, exp,
                 FALSE, debuglevel);

    }/*FOR*/

  }/*THEN*/

  UNPROTECT(1);

}/*C_PER_NODE_SCORE*/

