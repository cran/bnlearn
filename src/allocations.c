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

void *Calloc1D(size_t R, size_t size) {

void *p = NULL;

  if (R == 0)
    return NULL;

  p = calloc(R, size);

  if (!p)
    error("unable to allocate a %d array.", R);

  return(p);

}/*CALLOC1D*/

void BN_Free1D(void *p) {

  free(p);

}/*FREE1D*/

void **Calloc2D(size_t R, size_t C, size_t size) {

void **p = NULL;

  /* no corner cases, both dimensions required to be positive. */
  if ((R == 0) || (C == 0))
    error("trying to allocate a %dx%d two-dimensional array.", R, C);

  p = Calloc1D(R, sizeof(void *));

  for (int i = 0; i < R; i++)
    p[i] = Calloc1D(C, size);

  return p;

}/*CALLOC2D*/

void BN_Free2D(void **p, size_t R) {

int i = 0;

  for (i = 0; i < R; i++)
    free(p[i]);
  free(p);

}/*FREE2D*/

void ***Calloc3D(size_t R, size_t C, size_t L, size_t size) {

void ***p = NULL;

  /* no corner cases, all three dimensions required to be positive. */
  if ((R == 0) || (C == 0) || (L == 0))
    error("trying to allocate a %dx%dx%d three-dimensional array.", R, C, L);

  p = Calloc1D(R, sizeof(void *));
  for (int i = 0; i < R; i++)
    p[i] = Calloc2D(C, L, size);

  return p;

}/*CALLOC3D*/

void BN_Free3D(void ***p, size_t R, size_t C) {

int i = 0;

  for (i = 0; i < R; i++)
    BN_Free2D(p[i], C);
  free(p);

}/*FREE3D*/

