#include "../include/rcore.h"
#include "allocations.h"

void *Calloc1D(size_t R, size_t size) {

void *p = NULL;

  if (R == 0)
    return NULL;

  p = calloc(R, size);

  if (!p)
    error("unable to allocate a %llu array.", (unsigned long long)R);

  return p;

}/*CALLOC1D*/

void *Realloc1D(void *p, size_t R, size_t size) {

  p = realloc(p, R * size);

  if (!p)
    error("unable to reallocate a %llu array.", (unsigned long long)R);

  return p;

}/*REALLOC1D*/

void BN_Free1D(void *p) {

  free(p);

}/*FREE1D*/

void **Calloc2D(size_t R, size_t C, size_t size) {

void **p = NULL;

  /* no corner cases, both dimensions required to be positive. */
  if ((R == 0) || (C == 0))
    error("trying to allocate a %llux%llu two-dimensional array.",
        (unsigned long long)R, (unsigned long long)C);

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
    error("trying to allocate a %llux%llux%llu three-dimensional array.",
        (unsigned long long)R, (unsigned long long)C, (unsigned long long)L);

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

