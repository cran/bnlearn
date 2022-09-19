#include "../include/rcore.h"
#include "scores.h"

#define ENTRY(key, value) if (strcmp(label, key) == 0) return value;

score_e score_to_enum(const char *label) {

  ENTRY("loglik", LOGLIK);
  ENTRY("aic", AIC);
  ENTRY("bic", BIC);
  ENTRY("ebic", EBIC);
  ENTRY("bde", BDE);
  ENTRY("bds", BDS);
  ENTRY("bdj", BDJ);
  ENTRY("k2", K2);
  ENTRY("mbde", MBDE);
  ENTRY("bdla", BDLA);
  ENTRY("pred-loglik", PRED_LOGLIK);
  ENTRY("fnml", FNML);
  ENTRY("qnml", QNML);
  ENTRY("loglik-g", LOGLIK_G);
  ENTRY("aic-g", AIC_G);
  ENTRY("bic-g", BIC_G);
  ENTRY("ebic-g", EBIC_G);
  ENTRY("bge", BGE);
  ENTRY("pred-loglik-g", PRED_LOGLIK_G);
  ENTRY("loglik-cg", LOGLIK_CG);
  ENTRY("aic-cg", AIC_CG);
  ENTRY("bic-cg", BIC_CG);
  ENTRY("ebic-cg", EBIC_CG);
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

