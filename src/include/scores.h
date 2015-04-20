
/* score delta from score.delta.c */
SEXP score_delta(SEXP arc, SEXP network, SEXP data, SEXP score,
    SEXP score_delta, SEXP reference_score, SEXP op, SEXP extra, SEXP decomposable);

/* from graph.priors.c */
double graph_prior_prob(SEXP prior, SEXP target, SEXP cache, SEXP beta,
    int debuglevel);

/* from per.node.score.c */
SEXP per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, SEXP debug);
void c_per_node_score(SEXP network, SEXP data, SEXP score, SEXP targets,
    SEXP extra_args, int debuglevel, double *res);

/* score functions exported to per.node.score.c */
double dlik(SEXP x, double *nparams);
double cdlik(SEXP x, SEXP y, double *nparams);
double loglik_dnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel);
double glik(SEXP x, double *nparams);
double cglik(SEXP x, SEXP data, SEXP parents, double *nparams);
double c_fast_ccgloglik(double *xx, double **gp, int ngp, int nobs, int *config,
    int nconfig);
double loglik_gnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel);
double loglik_cgnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel);
double dirichlet_node(SEXP target, SEXP x, SEXP data, SEXP iss, SEXP prior,
    SEXP beta, SEXP experimental, int sparse, int debuglevel);
double wishart_node(SEXP target, SEXP x, SEXP data, SEXP isize, SEXP phi,
    SEXP prior, SEXP beta, int debuglevel);

