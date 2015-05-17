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

/* allocate a 2-dimensional contingency table. */
int **alloc2dcont (int length, int width) {

int **p = NULL, k = 0;

  p = (int **) R_alloc(length, sizeof(int *));

  for (k = 0; k < length; k++)
    p[k] = (int *) R_calloc(width, sizeof(int));

  return p;

}/*ALLOC2DCONT*/

/* allocate a 3-dimensional contingency table. */
int ***alloc3dcont (int length, int width, int depth) {

int ***p = NULL, i = 0, j = 0;

  p = (int ***) R_alloc(length, sizeof(int *));
  for (i = 0; i < length; i++) {

    p[i] = (int **) R_alloc(width, sizeof(int *));

    for (j = 0; j < width; j++)
      p[i][j] = (int *) R_calloc(depth, sizeof(int));

  }/*FOR*/

  return p;

}/*ALLOC3DCONT*/

/* allocate a 2-dimensional contingency table. */
double **alloc2dreal (int length, int width) {

double **p = NULL;
int k = 0;

  p = (double **) R_alloc(length, sizeof(double *));

  for (k = 0; k < length; k++)
    p[k] = (double *) R_calloc(width, sizeof(double));

  return p;

}/*ALLOC2DREAL*/

/* allocate a 3-dimensional real table. */
double ***alloc3dreal (int length, int width, int depth) {

double ***p = NULL;
int i = 0, j = 0;

  p = (double ***) R_alloc(length, sizeof(double *));
  for (i = 0; i < length; i++) {

    p[i] = (double **) R_alloc(width, sizeof(double *));

    for (j = 0; j < width; j++)
      p[i][j] = (double *) R_calloc(depth, sizeof(double));

  }/*FOR*/

  return p;

}/*ALLOC3DREAL*/

