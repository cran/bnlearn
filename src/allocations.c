#include <R.h>

/* allocate a 1-dimensional contingency table. */
int *alloc1dcont (int length) {

  int *p = (int *) R_alloc(length, sizeof(int));
  memset(p, '\0', sizeof(int) * length);

  return p;

}/*ALLOC1DCONT*/

/* allocate a 2-dimensional contingency table. */
int **alloc2dcont (int length, int width) {

  int **p;
  int k = 0;

  p = (int **) R_alloc(length, sizeof(int *));

  for (k = 0; k < length; k++) {

    p[k] = (int *) R_alloc(width, sizeof(int));
    memset(p[k], '\0', sizeof(int) * width);

  }/*FOR*/

  return p;

}/*ALLOC2DCONT*/

/* allocate a 3-dimensional contingency table. */
int ***alloc3dcont (int length, int width, int depth) {

  int ***p;
  int i = 0, j = 0;

  p = (int ***) R_alloc(length, sizeof(int *));
  for (i = 0; i < length; i++) {

    p[i] = (int **) R_alloc(width, sizeof(int *));

    for (j = 0; j < width; j++) {

      p[i][j] = (int*) R_alloc(depth, sizeof(int));
      memset(p[i][j], '\0', sizeof(int) * depth);

    }/*FOR*/

  }/*FOR*/

  return p;

}/*ALLOC3DCONT*/

/* allocate and initialize a status vector. */
short int *allocstatus (int length) {

  short int *p = (short int *) R_alloc(length, sizeof(short int));
  memset(p, '\0', sizeof(short int) * length);

  return p;

}/*ALLOCSTATUS*/

/* allocate a 1-dimensional real vector. */
double *alloc1dreal (int length) {

  double *p = (double *) R_alloc(length, sizeof(double));
  memset(p, '\0', sizeof(double) * length);

  return p;

}/*ALLOC1DREAL*/
