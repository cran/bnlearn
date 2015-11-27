#include "include/rcore.h"
#include "include/tests.h"

#define ENTRY(key, value) if (strcmp(label, key) == 0) return value;

test_e test_label(const char *label) {

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

}/*TEST_LABEL*/

