#include "include/rcore.h"
#include "include/scores.h"

#define DEBUG_BEFORE() \
  if (debugging) { \
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
    SEXP extra_args, bool debugging, double *res) {

int i = 0, ntargets = length(targets);
score_e s = score_to_enum(CHAR(STRING_ELT(score, 0)));
double nparams = 0, *k = NULL;
SEXP cur, iss, prior, beta, exp, l, nu, iss_w, newdata, custom_fn, custom_args;

  /* allocate dummy variable for the current node's label. */
  PROTECT(cur = allocVector(STRSXP, 1));

  switch(s) {

    /* discrete log-likelihood score. */
    case LOGLIK:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_dnode(cur, network, data, NULL, debugging);

      }/*FOR*/
      break;

    /* Gaussian log-likelihood score. */
    case LOGLIK_G:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_gnode(cur, network, data, NULL, debugging);

      }/*FOR*/
      break;

    /* Conditional Linear Gaussian log-likelihood score. */
    case LOGLIK_CG:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_cgnode(cur, network, data, NULL, debugging);

      }/*FOR*/
      break;

    /* AIC and BIC scores, discrete data. */
    case AIC:
    case BIC:

      k = REAL(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_dnode(cur, network, data, &nparams, debugging);
        res[i] -= (*k) * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

      }/*FOR*/
      break;

    /* AIC and BIC scores, Gaussian data. */
    case AIC_G:
    case BIC_G:

      k = REAL(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_gnode(cur, network, data, &nparams, debugging);
        res[i] -= (*k) * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

      }/*FOR*/
      break;

    /* AIC and BIC scores, Conditional Linear Gaussian data. */
    case AIC_CG:
    case BIC_CG:

      k = REAL(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_cgnode(cur, network, data, &nparams, debugging);
        res[i] -= (*k) * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", *k, nparams, (*k) * nparams);

      }/*FOR*/
      break;

    /* Bayesian Dirichlet equivalent score (BDe) and sparse score (BDs). */
    case BDE:
    case BDS:

      iss = getListElement(extra_args, "iss");
      prior = getListElement(extra_args, "prior");
      beta = getListElement(extra_args, "beta");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = dirichlet_node(cur, network, data, iss, FALSE, prior, beta,
                   R_NilValue, (s == BDS), debugging);

      }/*FOR*/
      break;

    /* Bayesian Dirichlet score with Jeffrey's prior, K2 score. */
    case BDJ:
    case K2:
      for (i = 0; i < ntargets; i++) {

        PROTECT(iss = (s == K2) ? ScalarReal(1) : ScalarReal(0.5));

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = dirichlet_node(cur, network, data, iss, TRUE,
                   R_NilValue, R_NilValue, R_NilValue, FALSE, debugging);

        UNPROTECT(1);

      }/*FOR*/
      break;

    /* Bayesian Gaussian equivalent score (BGe). */
    case BGE:

      iss = getListElement(extra_args, "iss.mu");
      nu = getListElement(extra_args, "nu");
      iss_w = getListElement(extra_args, "iss.w");
      prior = getListElement(extra_args, "prior");
      beta = getListElement(extra_args, "beta");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = wishart_node(cur, network, data, iss, nu, iss_w, prior,
                   beta, debugging);

      }/*FOR*/
      break;

    /* Mixture Bayesian Dirichlet equivalent score (mBDe). */
    case MBDE:

      iss = getListElement(extra_args, "iss");
      exp = getListElement(extra_args, "exp");
      prior = getListElement(extra_args, "prior");
      beta = getListElement(extra_args, "beta");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = dirichlet_node(cur, network, data, iss, FALSE, prior, beta,
                   exp, FALSE, debugging);

      }/*FOR*/
      break;

    /* Bayesian Dirichlet equivalent score, locally averaged. */
    case BDLA:

      prior = getListElement(extra_args, "prior");
      beta = getListElement(extra_args, "beta");
      l = getListElement(extra_args, "l");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = dirichlet_averaged_node(cur, network, data, l, prior,
                   beta, FALSE, debugging);

      }/*FOR*/
      break;

    /* discrete predictive log-likelihood. */
    case PRED_LOGLIK:

      newdata = getListElement(extra_args, "newdata");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = predictive_loglik_dnode(cur, network, data, newdata, &nparams,
                   debugging);

      }/*FOR*/
      break;

    /* Gaussian predictive log-likelihood. */
    case PRED_LOGLIK_G:

      newdata = getListElement(extra_args, "newdata");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = predictive_loglik_gnode(cur, network, data, newdata, &nparams,
                   debugging);

      }/*FOR*/
      break;

    case PRED_LOGLIK_CG:

      newdata = getListElement(extra_args, "newdata");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = predictive_loglik_cgnode(cur, network, data, newdata, &nparams,
                   debugging);

      }/*FOR*/
      break;

    case CUSTOM:

      custom_fn = getListElement(extra_args, "fun");
      custom_args = getListElement(extra_args, "args");

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = custom_score_function(cur, network, data, custom_fn,
                   custom_args, debugging);

      }/*FOR*/
      break;

    default:
      error("unknown score function.");

  }/*SWITCH*/

  UNPROTECT(1);

}/*C_PER_NODE_SCORE*/

