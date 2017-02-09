
/* enum for fitted node types, to be matched from the class in the R object. */
typedef enum {
  ENOFIT    =  0, /* error code, no such node type. */
  DNODE     =  1, /* categorical node */
  ONODE     =  2, /* ordinal node */
  GNODE     =  3, /* Gaussian node */
  CGNODE    =  4, /* conditional Gaussian node */
} fitted_node_e;

fitted_node_e r_fitted_node_label(SEXP class);
