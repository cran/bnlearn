#ifndef BN_FIT_OBJECTS_HEADER
#define BN_FIT_OBJECTS_HEADER

/* enum for fitted node types, to match the class in the R objects. */
typedef enum {
  ENOFIT    =  0, /* error code, no such node type. */
  DNODE     =  1, /* categorical node. */
  ONODE     =  2, /* ordinal node. */
  GNODE     =  3, /* Gaussian node. */
  CGNODE    =  4  /* conditional Gaussian node. */
} fitted_node_e;

/* enum for fitted network types, to match the class in the R objects. */
typedef enum {
  ENONET    =  0, /* error code, no such network type. */
  DNET      =  1, /* discrete Bayesian networks. */
  ONET      =  2, /* ordinal Bayesian networks. */
  DONET     =  3, /* mixed categorical and ordinal nodes. */
  GNET      =  4, /* Gaussian Bayesian networks. */
  CGNET     =  5  /* conditional Gaussian networks. */
} fitted_net_e;

fitted_node_e fitted_node_to_enum(SEXP object);
fitted_net_e fitted_net_to_enum(SEXP object);

/* a local distribution, meant to be embedded in the fitted_bn struct below. */
typedef struct {

  int nparents;      /* number of parents. */
  int *parents;      /* indexes of the parent nodes. */

  union {

    struct {

      int ndims;       /* number of dimensions of the CPT. */
      int *dims;       /* dimensions of the CPT. */
      int nconfigs;    /* number of parents configurations. */
      double *cpt;     /* conditional probability table. */

    } d;

    struct {

      int ncoefs;      /* number of regression coefficients. */
      double *coefs;   /* regression coefficients. */
      double sd;       /* standard error of the residuals. */

    } g;

    struct {

      int ndparents;   /* number of discrete parents. */
      int ngparents;   /* number of continuous parents. */
      int *dparents;   /* indexes of the discrete parents in the parents array. */
      int *gparents;   /* indexes of the continuous parents in the parents array. */
      int ncoefs;      /* number of regression coefficients. */
      int nconfigs;    /* number of discrete parents configurations. */
      double *coefs;   /* regression coefficients. */
      double *sd;      /* standard errors of the residuals. */

    } cg;

  };

} ldist;

/* a fitted Bayesian network, mapping to the bn.fit class in R code. */
typedef struct {

  fitted_net_e type;           /* network type (discrete, Gaussian, etc.). */
  int nnodes;                  /* number of nodes in the network. */
  const char **labels;         /* node labels. */
  fitted_node_e *node_types;   /* node types (discrete, Gaussian, etc.). */
  ldist *ldists;               /* local distributions for the nodes. */

} fitted_bn;

fitted_bn fitted_network_from_SEXP(SEXP fitted);
void print_fitted_network(fitted_bn);
void FreeFittedBN(fitted_bn bn);

double nparams_fitted_node(ldist ld, fitted_node_e type, bool effective);

#endif
