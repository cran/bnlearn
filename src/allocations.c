#include "include/rcore.h"

static int stack_counter = 0;

/* instrumentation to debug PROTECT()/UNPROTECT() calls. */
void PROTECT_DEBUG(SEXP s, const char *fun, const char *file, int line) {

  Rprintf("[%s()][%d -> %d] PROTECT() at %s:%d\n", fun,
    stack_counter, stack_counter + 1, file, line);
  stack_counter++;
  Rf_protect(s);

}/*PROTECT_DEBUG*/

void UNPROTECT_DEBUG(int n, const char *fun, const char *file, int line) {

  Rprintf("[%s()][%d -> %d] UNPROTECT() at %s:%d\n", fun,
    stack_counter, stack_counter - n, file, line);
  stack_counter -= n;
  Rf_unprotect(n);

}/*UNPROTECT_DEBUG*/

/* make R_alloc() behave like calloc() and zero allocated memory. */
void *R_calloc(size_t num, size_t size) {

void *p = R_alloc(num, size);

  /* memset()ting a NULL pointer is a NOP. */
  if (p)
    memset(p, '\0', size * num);
  else if (num > 0)
    error("memsetting() %d bytes on a NULL pointer.", size * num);

  return p;

}/*R_CALLOC*/

/* two-dimensional version of R_calloc(). */
void **R_calloc2(size_t num1, size_t num2, size_t size) {

void **p = NULL;

  /* no corner cases, both dimensions required to be positive. */
  if ((num1 == 0) || (num2 == 0))
    error("trying to allocate a %dx%d two-dimensional array.", num1, num2);

  p = (void **) R_alloc(num1, sizeof(void *));

  for (int i = 0; i < num1; i++)
    p[i] = R_calloc(num2, size);

  return p;

}/*R_CALLOC2*/

/* three-dimensional version of R_calloc(). */
void ***R_calloc3(size_t num1, size_t num2, size_t num3, size_t size) {

void ***p = NULL;

  /* no corner cases, all three dimensions required to be positive. */
  if ((num1 == 0) || (num2 == 0) || (num3 == 0))
    error("trying to allocate a %dx%dx%d three-dimensional array.",
      num1, num2, num3);

  p = (void ***) R_alloc(num1, sizeof(void *));
  for (int i = 0; i < num1; i++)
    p[i] = R_calloc2(num2, num3, size);

  return p;

}/*R_CALLOC3*/

