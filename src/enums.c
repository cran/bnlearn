#include "include/rcore.h"
#include "include/tests.h"
#include "include/scores.h"
#include "include/fitted.h"
#include "include/preprocessing.h"

#define ENTRY(key, value) if (strcmp(label, key) == 0) return value;

test_e test_to_enum(const char *label) {

  ENTRY("mi", MI);
  ENTRY("mi-adf", MI_ADF);
  ENTRY("x2", X2);
  ENTRY("x2-adf", X2_ADF);
  ENTRY("jt", JT);
  ENTRY("cor", COR);
  ENTRY("zf", ZF);
  ENTRY("mi-g", MI_G);
  ENTRY("mi-cg", MI_CG);
  ENTRY("mi-sh", MI_SH);
  ENTRY("mi-g-sh", MI_G_SH);
  ENTRY("mc-mi", MC_MI);
  ENTRY("mc-x2", MC_X2);
  ENTRY("sp-mi", SP_MI);
  ENTRY("sp-x2", SP_X2);
  ENTRY("mc-jt", MC_JT);
  ENTRY("smc-mi", SMC_MI);
  ENTRY("smc-x2", SMC_X2);
  ENTRY("smc-jt", SMC_JT);
  ENTRY("mc-cor", MC_COR);
  ENTRY("mc-mi-g", MC_MI_G);
  ENTRY("mc-zf", MC_ZF);
  ENTRY("smc-cor", SMC_COR);
  ENTRY("smc-zf", SMC_ZF);
  ENTRY("smc-mi-g", SMC_MI_G);

  return ENOTEST;

}/*TEST_TO_ENUM*/

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

score_e score_to_enum(const char *label) {

  ENTRY("loglik", LOGLIK);
  ENTRY("aic", AIC);
  ENTRY("bic", BIC);
  ENTRY("bde", BDE);
  ENTRY("bds", BDS);
  ENTRY("bdj", BDJ);
  ENTRY("k2", K2);
  ENTRY("mbde", MBDE);
  ENTRY("bdla", BDLA);
  ENTRY("pred-loglik", PRED_LOGLIK);
  ENTRY("loglik-g", LOGLIK_G);
  ENTRY("aic-g", AIC_G);
  ENTRY("bic-g", BIC_G);
  ENTRY("bge", BGE);
  ENTRY("pred-loglik-g", PRED_LOGLIK_G);
  ENTRY("loglik-cg", LOGLIK_CG);
  ENTRY("aic-cg", AIC_CG);
  ENTRY("bic-cg", BIC_CG);
  ENTRY("pred-loglik-cg", PRED_LOGLIK_CG);
  ENTRY("custom", CUSTOM);

  return ENOSCORE;

}/*SCORE_TO_ENUM*/

gprior_e gprior_to_enum(const char *label) {

  ENTRY("uniform", UNIFORM);
  ENTRY("vsp", VSP);
  ENTRY("cs", CS);
  ENTRY("marginal", MU);

  return ENOPRIOR;

}/*GPRIOR_TO_ENUM*/

discretization_e discretization_to_enum(const char *label) {

  ENTRY("quantile", QUANTILE);
  ENTRY("interval", INTERVAL);
  ENTRY("hartemink", HARTEMINK);

  return ENOMETHOD;

}/*DISCRETIZATION_TO_ENUM*/
