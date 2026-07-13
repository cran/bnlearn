#include "../include/rcore.h"
#include "parameters.h"

#define ENTRY(key, value) \
  do { \
    if (strcmp(label, key) == 0) return value; \
  } while (0)

estimator_e estimator_to_enum(const char *label) {

  ENTRY("mle", MLE);
  ENTRY("bayes", BAYES);
  ENTRY("hdir", HDIR);
  ENTRY("hard-em", HARD_EM);
  ENTRY("mle-g", MLE_G);
  ENTRY("hard-em-g", HARD_EM_G);
  ENTRY("mle-cg", MLE_CG);
  ENTRY("hard-em-cg", HARD_EM_CG);
  ENTRY("mle-zihp", MLE_ZIHP);
  ENTRY("mle-zinb", MLE_ZINB);

  return ENOEST;

}/*ESTIMATOR_TO_ENUM*/
