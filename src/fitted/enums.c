#include "../include/rcore.h"
#include "../minimal/common.h"
#include "fitted.h"

fitted_node_e fitted_node_to_enum(SEXP object) {

  if (c_is(object, "bn.fit.dnode"))
    return DNODE;
  else if (c_is(object, "bn.fit.onode"))
    return ONODE;
  else if (c_is(object, "bn.fit.gnode"))
    return GNODE;
  else if (c_is(object, "bn.fit.cgnode"))
    return CGNODE;

  return ENOFIT;

}/*FITTED_NODE_TO_ENUM*/

fitted_net_e fitted_net_to_enum(SEXP object) {

  if (c_is(object, "bn.fit.dnet"))
    return DNET;
  else if (c_is(object, "bn.fit.onet"))
    return ONET;
  else if (c_is(object, "bn.fit.donet"))
    return DONET;
  else if (c_is(object, "bn.fit.gnet"))
    return GNET;
  else if (c_is(object, "bn.fit.cgnet"))
    return CGNET;

  return ENONET;

}/*FITTED_NET_TO_ENUM*/
