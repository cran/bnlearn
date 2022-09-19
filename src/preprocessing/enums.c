#include "../include/rcore.h"
#include "preprocessing.h"

#define ENTRY(key, value) if (strcmp(label, key) == 0) return value;

discretization_e discretization_to_enum(const char *label) {

  ENTRY("quantile", QUANTILE);
  ENTRY("interval", INTERVAL);
  ENTRY("hartemink", HARTEMINK);

  return ENOMETHOD;

}/*DISCRETIZATION_TO_ENUM*/
