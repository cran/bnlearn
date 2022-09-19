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

fitted_node_e fitted_node_to_enum(SEXP class);
fitted_net_e fitted_net_to_enum(SEXP class);

#endif
