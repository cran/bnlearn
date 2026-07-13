#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/common.h"
#include "scores.h"

#define DEBUG_BEFORE() \
  do { \
  if (debugging) { \
    Rprintf("----------------------------------------------------------------\n"); \
    Rprintf("* processing node %s.\n", CHAR(STRING_ELT(cur, 0))); \
  }/*THEN*/ \
  } while (0)

/* read the EM controls (maximum iterations, convergence tolerance, M-step type)
 * for the EM count scores from the extra arguments, falling back to the same
 * defaults as bn.fit() when an argument is absent. */
static void get_em_controls(SEXP extra_args, int *em_max_iter, double *em_tol,
    bool *one_step) {

SEXP mi = getListElement(extra_args, "em.max.iter");
SEXP tol = getListElement(extra_args, "em.tol");
SEXP ms = getListElement(extra_args, "m.step");

  *em_max_iter = isNull(mi) ? 100 : asInteger(mi);
  *em_tol = isNull(tol) ? 1e-8 : asReal(tol);
  *one_step = !isNull(ms) && (strcmp(CHAR(STRING_ELT(ms, 0)), "one-step") == 0);

}/*GET_EM_CONTROLS*/

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

int i = 0, ntargets = length(targets), nparents = 0;
score_e s = score_to_enum(CHAR(STRING_ELT(score, 0)));
double nparams = 0, k = 0, *gamma = NULL;
SEXP cur, iss, prior, beta, exp, nu, iss_w, newdata, custom_fn, custom_args;

  /* allocate dummy variable for the current node's label. */
  PROTECT(cur = allocVector(STRSXP, 1));

  switch(s) {

    /* discrete log-likelihood score. */
    case LOGLIK:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_dnode(cur, network, data, NULL, NULL, debugging);

      }/*FOR*/
      break;

    /* Gaussian log-likelihood score. */
    case LOGLIK_G:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_gnode(cur, network, data, NULL, NULL, debugging);

      }/*FOR*/
      break;

    /* Conditional Linear Gaussian log-likelihood score. */
    case LOGLIK_CG:
      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_cgnode(cur, network, data, NULL, NULL, debugging);

      }/*FOR*/
      break;

    /* AIC and BIC scores, discrete data. */
    case AIC:
    case BIC:

      k = NUM(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_dnode(cur, network, data, &nparams, NULL, debugging);
        res[i] -= k * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/
      break;

    /* extended BIC score, discrete data. */
    case EBIC:

      k = NUM(getListElement(extra_args, "k"));
      gamma = REAL(getListElement(extra_args, "gamma"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_dnode(cur, network, data, &nparams, &nparents, debugging);
        res[i] -= k * nparams +
                    4 * (*gamma) * (nparents) * log((double)(length(data)));

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/
      break;

    /* AIC and BIC scores, Gaussian data. */
    case AIC_G:
    case BIC_G:

      k = NUM(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_gnode(cur, network, data, &nparams, NULL, debugging);
        res[i] -= k * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/
      break;

    /* extended BIC score, Gaussian data. */
    case EBIC_G:

      k = NUM(getListElement(extra_args, "k"));
      gamma = REAL(getListElement(extra_args, "gamma"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_gnode(cur, network, data, &nparams, &nparents, debugging);
        res[i] -= k * nparams +
                    4 * (*gamma) * (nparents) * log((double)(length(data)));

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/
      break;

    /* AIC and BIC scores, Conditional Linear Gaussian data. */
    case AIC_CG:
    case BIC_CG:

      k = NUM(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_cgnode(cur, network, data, &nparams, NULL, debugging);
        res[i] -= k * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/
      break;

    /* extended BIC score, Conditional Linear Gaussian data. */
    case EBIC_CG:

      k = NUM(getListElement(extra_args, "k"));
      gamma = REAL(getListElement(extra_args, "gamma"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = loglik_cgnode(cur, network, data, &nparams, &nparents, debugging);
        res[i] -= k * nparams +
                    4 * (*gamma) * (nparents) * log((double)(length(data)));

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

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

    /* Factorized Normalized Maximum Likelihood. */
    case FNML:

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = fnml_node(cur, network, data, debugging);

      }/*FOR*/
      break;

    /* Quotient Normalized Maximum Likelihood. */
    case QNML:

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = qnml_node(cur, network, data, debugging);

      }/*FOR*/

      Free1D(regret_table);

      break;

    /* Discrete (penalized) node-average likelihood. */
    case NAL:
    case PNAL:

      if (s == PNAL)
        k = NUM(getListElement(extra_args, "k"));
      else
        k = 0;

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = nal_dnode(cur, network, data, k, debugging);

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/

      break;

    /* Gaussian (penalized) node-average likelihood. */
    case NAL_G:
    case PNAL_G:

      if (s == PNAL_G)
        k = NUM(getListElement(extra_args, "k"));
      else
        k = 0;

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = nal_gnode(cur, network, data, k, debugging);

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/

      break;

    /* Gaussian (penalized) node-average likelihood. */
    case NAL_CG:
    case PNAL_CG:

      if (s == PNAL_CG)
        k = NUM(getListElement(extra_args, "k"));
      else
        k = 0;

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = nal_cgnode(cur, network, data, k, debugging);

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/

      break;

    /* zero-inflated hyper-Poisson / negative binomial likelihood, estimated by
     * EM with the GLM components fitted by IRLS. */
    case LOGLIK_ZIHP:
    case LOGLIK_ZINB: {

      int em_max_iter = 0;
      double em_tol = 0;
      bool one_step = FALSE;
      glm_family_e family = (s == LOGLIK_ZIHP) ? GLM_HYPERPOISSON : GLM_NEGBIN;
      get_em_controls(extra_args, &em_max_iter, &em_tol, &one_step);

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = em_irls_node(cur, network, data, NULL, NULL, FALSE, 0,
                   debugging, em_max_iter, em_tol, one_step, family);

      }/*FOR*/

      break;

    }/*CASE*/

    /* zero-inflated hyper-Poisson / negative binomial penalised likelihood. */
    case AIC_ZIHP:
    case BIC_ZIHP:
    case AIC_ZINB:
    case BIC_ZINB: {

      int em_max_iter = 0;
      double em_tol = 0;
      bool one_step = FALSE;
      glm_family_e family = (s == AIC_ZIHP || s == BIC_ZIHP) ? GLM_HYPERPOISSON : GLM_NEGBIN;
      get_em_controls(extra_args, &em_max_iter, &em_tol, &one_step);
      k = NUM(getListElement(extra_args, "k"));

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = em_irls_node(cur, network, data, NULL, &nparams, FALSE, 0,
                   debugging, em_max_iter, em_tol, one_step, family);
        res[i] -= k * nparams;

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/

      break;

    }/*CASE*/

    /* zero-inflated hyper-Poisson / negative binomial (node average) likelihood,
     * penalised or not. */
    case NAL_ZIHP:
    case PNAL_ZIHP:
    case NAL_ZINB:
    case PNAL_ZINB: {

      int em_max_iter = 0;
      double em_tol = 0;
      bool one_step = FALSE;
      glm_family_e family = (s == NAL_ZIHP || s == PNAL_ZIHP) ? GLM_HYPERPOISSON : GLM_NEGBIN;
      get_em_controls(extra_args, &em_max_iter, &em_tol, &one_step);

      if (s == PNAL_ZIHP || s == PNAL_ZINB)
        k = NUM(getListElement(extra_args, "k"));
      else
        k = 0;

      for (i = 0; i < ntargets; i++) {

        SET_STRING_ELT(cur, 0, STRING_ELT(targets, i));
        DEBUG_BEFORE();
        res[i] = em_irls_node(cur, network, data, NULL, &nparams, TRUE, k,
                   debugging, em_max_iter, em_tol, one_step, family);

        if (debugging)
          Rprintf("  > penalty is %lf x %.0lf = %lf.\n", k, nparams, k * nparams);

      }/*FOR*/

      break;

    }/*CASE*/

    /* custom, user-defined score. */
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

