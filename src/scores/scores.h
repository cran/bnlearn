#ifndef NETWORK_SCORES_HEADER
#define NETWORK_SCORES_HEADER

#include "../core/contingency.tables.h"

/* enum for scores, to be matched from the label string passed down from R. */
typedef enum {
  ENOSCORE       =   0, /* error code, no such score. */

  /* scores for discrete data. */
  LOGLIK         =   1, /* log-likelihood. */
  PRED_LOGLIK    =   2, /* predictive log-likelihood. */
  AIC            =   3, /* AIC. */
  BIC            =   4, /* BIC. */
  EBIC           =   5, /* extended BIC. */
  BDE            =   6, /* Bayesian Dirichlet equivalent score. */
  BDS            =   7, /* Bayesian Dirichlet sparse score. */
  BDJ            =   8, /* Bayesian Dirichlet with Jeffrey's prior. */
  K2             =   9, /* K2 score. */
  MBDE           =  10, /* Bayesian Dirichlet equivalent score, interventional data .*/
  BDLA           =  11, /* Bayesian Dirichlet score, locally averaged. */
  FNML           =  12, /* Factorized Normalized Maximum Likelihood. */
  QNML           =  13, /* Quotient Normalized Maximum Likelihood. */
  NAL            =  14, /* Node-average likelihood, discrete data. */
  PNAL           =  15, /* Penalized node-average likelihood. */

  /* scores for Gaussian data. */
  LOGLIK_G       = 100, /* log-likelihood. */
  PRED_LOGLIK_G  = 101, /* predictive log-likelihood. */
  AIC_G          = 102, /* AIC. */
  BIC_G          = 103, /* BIC. */
  EBIC_G         = 104, /* extended BIC. */
  BGE            = 105, /* Bayesian Gaussian equivalent score. */
  NAL_G          = 106, /* Node-average likelihood. */
  PNAL_G         = 107, /* Penalized node-average likelihood. */

  /* scores for conditional Gaussian data. */
  LOGLIK_CG      = 200, /* log-likelihood. */
  PRED_LOGLIK_CG = 201, /* predictive log-likelihood. */
  AIC_CG         = 202, /* AIC. */
  BIC_CG         = 203, /* BIC. */
  EBIC_CG        = 204, /* extended BIC. */
  NAL_CG         = 205, /* Node-average likelihood. */
  PNAL_CG        = 206, /* Penalized node-average likelihood. */

  CUSTOM         = 300  /* user-provided score function. */
} score_e;

score_e score_to_enum(const char *label);

/* enum for graph priors, to be matched from the label string passed down from R. */
typedef enum {
  ENOPRIOR  =  0, /* error code, no such graph prior. */
  UNIFORM   =  1, /* uniform prior. */
  VSP       =  2, /* variable selection prior. */
  CS        =  3, /* Castelo & Siebes prior. */
  MU        =  4, /* marginal uniform prior. */
} gprior_e;

gprior_e gprior_to_enum(const char *label);

/* score delta from score.delta.c */
SEXP score_delta(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra, SEXP decomposable);

/* from graph.priors.c */
double graph_prior_prob(SEXP prior, SEXP target, SEXP beta, SEXP cache,
    bool debugging);

/* from per.node.score.c */
SEXP per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, SEXP debug);
void c_per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, bool debugging, double *res);

/* score functions exported to per.node.score.c and to other score functions */
double dlik(counts1d marginal);
double cdlik(counts2d joint);
double loglik_dnode_root(SEXP x, double *nparams);
double loglik_dnode_parents(SEXP x, SEXP y, double *nparams);
double loglik_dnode(SEXP target, SEXP x, SEXP data, double *nparams,
    int *nparents, bool debugging);
double pdnode(SEXP x, SEXP new_x, double *nparams);
double cpdnode(SEXP x, SEXP y, SEXP x2, SEXP y2, double *nparams);
double predictive_loglik_dnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, bool debugging);
double glik(SEXP x, double *nparams);
double cglik(SEXP x, SEXP data, SEXP parents, double *nparams);
double pgnode(SEXP x, SEXP new_x, double *nparams);
double cpgnode(SEXP x, SEXP x2, SEXP data, SEXP newdata, SEXP parents,
    double *nparams);
double predictive_loglik_gnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, bool debugging);
double c_fast_ccgloglik(double *xx, double **gp, int ngp, int nobs, int *config,
    int nconfig);
double loglik_gnode(SEXP target, SEXP x, SEXP data, double *nparams,
    int *nparents, bool debugging);
double loglik_cgnode(SEXP target, SEXP x, SEXP data, double *nparams,
    int *np, bool debugging);
double predictive_loglik_cgnode(SEXP target, SEXP x, SEXP data, SEXP newdata,
    double *nparams, bool debugging);
double dirichlet_node(SEXP target, SEXP x, SEXP data, SEXP iss, int per_node,
    SEXP prior, SEXP beta, SEXP experimental, int sparse, bool debugging);
double dirichlet_averaged_node(SEXP target, SEXP x, SEXP data, SEXP l,
    SEXP prior, SEXP beta, int sparse, bool debugging);
double wishart_node(SEXP target, SEXP x, SEXP data, SEXP iss, SEXP nu,
    SEXP iss_w, SEXP prior, SEXP beta, bool debugging);
double custom_score_function(SEXP target, SEXP x, SEXP data, SEXP custom_fn,
    SEXP custom_args, bool debugging);
double fnml_node(SEXP target, SEXP x, SEXP data, bool debugging);
double qnml_node(SEXP target, SEXP x, SEXP data, bool debugging);
double nal_dnode_root(SEXP x, double k);
double nal_dnode_parents(SEXP x, SEXP y, double k);
double nal_dnode(SEXP target, SEXP x, SEXP data, double k, bool debugging);
double glik_incomplete(SEXP x, double k);
double cglik_incomplete(SEXP x, SEXP data, SEXP parents, double k);
double nal_gnode(SEXP target, SEXP x, SEXP data, double k, bool debugging);
double nal_cgnode(SEXP target, SEXP x, SEXP data, double k, bool debugging);

/* exports for the regret table of normalized maximum likelihood scores. */
#define MAX_REGRET_TABLE_N 1000
#define MAX_REGRET_TABLE_K 100
extern double *regret_table;

double nml_regret(double n, double k);
double *get_regret_table(int N, int K);

#endif
