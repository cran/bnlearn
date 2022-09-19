#include "../include/rcore.h"
#include "../minimal/common.h"
#include "fitted.h"

fitted_node_e fitted_node_to_enum(SEXP class) {

  if (c_is(class, "bn.fit.dnode"))
    return DNODE;
  else if (c_is(class, "bn.fit.onode"))
    return ONODE;
  else if (c_is(class, "bn.fit.gnode"))
    return GNODE;
  else if (c_is(class, "bn.fit.cgnode"))
    return CGNODE;

  return ENOFIT;

}/*FITTED_NODE_TO_ENUM*/

fitted_net_e fitted_net_to_enum(SEXP class) {

  if (c_is(class, "bn.fit.dnet"))
    return DNET;
  else if (c_is(class, "bn.fit.onet"))
    return ONET;
  else if (c_is(class, "bn.fit.donet"))
    return DONET;
  else if (c_is(class, "bn.fit.gnet"))
    return GNET;
  else if (c_is(class, "bn.fit.cgnet"))
    return CGNET;

  return ENONET;

}/*FITTED_NET_TO_ENUM*/
