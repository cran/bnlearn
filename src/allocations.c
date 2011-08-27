#include <R.h>

/* allocate a 1-dimensional contingency table. */
int *alloc1dcont (int length) {

int *p = (int *) R_alloc(length, sizeof(int));
memset(p, '\0', sizeof(int) * length);

  return p;

}/*ALLOC1DCONT*/

/* allocate a 2-dimensional contingency table. */
int **alloc2dcont (int length, int width) {

int **p = NULL, k = 0;

  p = (int **) R_alloc(length, sizeof(int *));

  for (k = 0; k < length; k++) {

    p[k] = (int *) R_alloc(width, sizeof(int));
    memset(p[k], '\0', sizeof(int) * width);

  }/*FOR*/

  return p;

}/*ALLOC2DCONT*/

/* allocate a 3-dimensional contingency table. */
int ***alloc3dcont (int length, int width, int depth) {

int ***p = NULL, i = 0, j = 0;

  p = (int ***) R_alloc(length, sizeof(int *));
  for (i = 0; i < length; i++) {

    p[i] = (int **) R_alloc(width, sizeof(int *));

    for (j = 0; j < width; j++) {

      p[i][j] = (int *) R_alloc(depth, sizeof(int));
      memset(p[i][j], '\0', sizeof(int) * depth);

    }/*FOR*/

  }/*FOR*/

  return p;

}/*ALLOC3DCONT*/

/* allocate and initialize a status vector. */
short int *allocstatus (int length) {

short int *p = NULL;

  p = (short int *) R_alloc(length, sizeof(short int));
  memset(p, '\0', sizeof(short int) * length);

  return p;

}/*ALLOCSTATUS*/

/* allocate a 1-dimensional real vector. */
double *alloc1dreal (int length) {

double *p = NULL;

  p = (double *) R_alloc(length, sizeof(double));
  memset(p, '\0', sizeof(double) * length);

  return p;

}/*ALLOC1DREAL*/

/* allocate a 2-dimensional contingency table. */
double **alloc2dreal (int length, int width) {

double **p = NULL;
int k = 0;

  p = (double **) R_alloc(length, sizeof(double *));

  for (k = 0; k < length; k++) {

    p[k] = (double *) R_alloc(width, sizeof(double));
    memset(p[k], '\0', sizeof(double) * width);

  }/*FOR*/

  return p;

}/*ALLOC2DREAL*/

/* allocate a 3-dimensional real table. */
double ***alloc3dreal (int length, int width, int depth) {

double ***p = NULL;
int i = 0, j = 0;

  p = (double ***) R_alloc(length, sizeof(double *));
  for (i = 0; i < length; i++) {

    p[i] = (double **) R_alloc(width, sizeof(double *));

    for (j = 0; j < width; j++) {

      p[i][j] = (double *) R_alloc(depth, sizeof(double));
      memset(p[i][j], '\0', sizeof(double) * depth);

    }/*FOR*/

  }/*FOR*/

  return p;

}/*ALLOC3DREAL*/

void **alloc1dpointer (int length) {

void **p = NULL;

  p = (void **) R_alloc(length, sizeof(void *));

  return p;

}/*ALLOC1DPOINTERS*/

#define alloc1dstring(length) ((char **) alloc1dpointer(length))

