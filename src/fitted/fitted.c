#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../minimal/common.h"
#include "fitted.h"

/* create a fitted_bn data structure from a SEXP, copying as little as possible. */
fitted_bn fitted_network_from_SEXP(SEXP fitted) {

fitted_net_e type = fitted_net_to_enum(fitted);
fitted_bn bn;
SEXP temp, nodes_in_fitted_bn, cur_node, parents;

  /* find out the type and the size of the network, and allocate its
   * components (node types and labels, local distributions. */
  bn.type = type;
  bn.nnodes = length(fitted);
  bn.node_types = Calloc1D(bn.nnodes, sizeof(fitted_node_e));
  bn.ldists = Calloc1D(bn.nnodes, sizeof(ldist));
  bn.labels = Calloc1D(bn.nnodes, sizeof(const char *));

  /* copy the node labels. */
  PROTECT(nodes_in_fitted_bn = getAttrib(fitted, R_NamesSymbol));
  for (int i = 0; i < bn.nnodes; i++)
    bn.labels[i] = CHAR(STRING_ELT(nodes_in_fitted_bn, i));

  /* save pointers to the parameters in the SEXP and their dimensions. */
  for (int i = 0; i < bn.nnodes; i++) {

    cur_node = VECTOR_ELT(fitted, i);
    bn.node_types[i] = fitted_node_to_enum(cur_node);

    temp = getListElement(cur_node, "parents");
    PROTECT(parents = match(nodes_in_fitted_bn, temp, 0));
    bn.ldists[i].nparents = length(parents);
    bn.ldists[i].parents = Calloc1D(length(parents), sizeof(int));
    for (int j = 0; j < length(parents); j++)
      bn.ldists[i].parents[j] = INTEGER(parents)[j] - 1;
    UNPROTECT(1);

    if ((bn.node_types[i] == DNODE) || (bn.node_types[i] == ONODE)) {

      /* discrete nodes. */
      temp = getListElement(cur_node, "prob");
      bn.ldists[i].d.cpt = REAL(temp);
      bn.ldists[i].d.nconfigs = length(temp);
      temp = getAttrib(temp, R_DimSymbol);
      bn.ldists[i].d.ndims = length(temp);
      bn.ldists[i].d.dims = INTEGER(temp);
      /* for root nodes, the number of parent configurations will be 1. */
      bn.ldists[i].d.nconfigs /= bn.ldists[i].d.dims[0];

    }/*THEN*/
    else if (bn.node_types[i] == GNODE) {

      /* Gaussian nodes. */
      temp = getListElement(cur_node, "coefficients");
      bn.ldists[i].g.ncoefs = length(temp);
      bn.ldists[i].g.coefs = REAL(temp);
      temp = getListElement(cur_node, "sd");
      bn.ldists[i].g.sd = NUM(temp);

    }/*THEN*/
    else if (bn.node_types[i] == CGNODE) {

      /* conditional Gaussian nodes. */
      temp = getListElement(cur_node, "dparents");
      bn.ldists[i].cg.ndparents = length(temp);
      bn.ldists[i].cg.dparents = Calloc1D(length(temp), sizeof(int));
      for (int j = 0; j < length(temp); j++)
        bn.ldists[i].cg.dparents[j] = bn.ldists[i].parents[INTEGER(temp)[j] - 1];
      temp = getListElement(cur_node, "gparents");
      bn.ldists[i].cg.ngparents = length(temp);
      bn.ldists[i].cg.gparents = Calloc1D(length(temp), sizeof(int));
      for (int j = 0; j < length(temp); j++)
        bn.ldists[i].cg.gparents[j] = bn.ldists[i].parents[INTEGER(temp)[j] - 1];
      temp = getListElement(cur_node, "coefficients");
      bn.ldists[i].cg.ncoefs = INTEGER(getAttrib(temp, R_DimSymbol))[0];
      bn.ldists[i].cg.nconfigs = INTEGER(getAttrib(temp, R_DimSymbol))[1];
      bn.ldists[i].cg.coefs = REAL(temp);
      temp = getListElement(cur_node, "sd");
      bn.ldists[i].cg.sd = REAL(temp);

    }/*THEN*/

  }/*FOR*/

  UNPROTECT(1);

  return bn;

}/*FITTED_NETWORK_FROM_SEXP*/

/* print a short summary of a fitted_bn data structure. */
void print_fitted_network(fitted_bn bn) {

int parconfigs = 0;

  switch(bn.type) {

    case DNET:
    case ONET:
    case DONET:
      Rprintf("discrete network: ");
      break;

    case GNET:
      Rprintf("Gaussian network: ");
      break;

    case CGNET:
      Rprintf("conditional Gaussian network: ");
      break;

    case ENONET:
    default:
      Rprintf("unknown network type: ");

  }/*SWTICH*/

  Rprintf("%d nodes.\n", bn.nnodes);

  for (int i = 0; i < bn.nnodes; i++) {

    Rprintf("%10s", bn.labels[i]);

    switch(bn.node_types[i]) {

      case DNODE:
      case ONODE:

        Rprintf(" [%s ]: %d parents, %d levels",
            (bn.node_types[i] == DNODE) ? "D" : "O",
            bn.ldists[i].nparents, bn.ldists[i].d.dims[0]);

        parconfigs = 1;
        if (bn.ldists[i].nparents > 0) {

          for (int j = 1; j < bn.ldists[i].d.ndims; j++)
            parconfigs *= bn.ldists[i].d.dims[j];

          Rprintf(", %d configurations.\n", parconfigs);

        }/*THEN*/
        else {

          Rprintf(".\n");

        }/*ELSE*/

        Rprintf("         CPT: ");
        for (int j = 0; j < MIN(bn.ldists[i].d.dims[0] * parconfigs, 5); j++)
          Rprintf("%g ", bn.ldists[i].d.cpt[j]);
        Rprintf("\n");

        break;

      case GNODE:
        Rprintf(" [G ]: %d parents.\n", bn.ldists[i].nparents);
        Rprintf("         COEFS: ");
        for (int j = 0; j < MIN(bn.ldists[i].g.ncoefs, 5); j++)
          Rprintf("%g ", bn.ldists[i].g.coefs[j]);
        Rprintf("SD: %g\n", bn.ldists[i].g.sd);
        break;

      case CGNODE:
        Rprintf(" [CG]: %d parents, %d coefficients, %d configurations.\n",
            bn.ldists[i].nparents, bn.ldists[i].cg.ncoefs,
            bn.ldists[i].cg.nconfigs);
        Rprintf("         COEFS: ");
        for (int j = 0; j < MIN(bn.ldists[i].cg.ncoefs, 5); j++)
          Rprintf("%g ", bn.ldists[i].cg.coefs[j]);
        Rprintf("\n         SD: ");
        for (int j = 0; j < MIN(bn.ldists[i].cg.nconfigs, 5); j++)
          Rprintf("%g ", bn.ldists[i].cg.sd[j]);
        break;

      case ENOFIT:
      default:
        break;

    }/*SWTICH*/

  }/*FOR*/

}/*PRINT_FITTED_NETWORK*/

/* free a fitted_bn data structure, including the local distributions. */
void FreeFittedBN(fitted_bn bn) {

  for (int i = 0; i < bn.nnodes; i++)
    Free1D(bn.ldists[i].parents);
  for (int i = 0; i < bn.nnodes; i++) {

    if (bn.node_types[i] == CGNODE) {

      Free1D(bn.ldists[i].cg.dparents);
      Free1D(bn.ldists[i].cg.gparents);

    }/*THEN*/

  }/*THEN*/
  Free1D(bn.node_types);
  Free1D(bn.ldists);
  Free1D(bn.labels);

}/*FREEFITTEDBN*/

