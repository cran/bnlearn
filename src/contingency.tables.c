#include "include/rcore.h"

/* initialize a one-dimensional contingency table. */
void fill_1d_table(int *xx, int **n, int llx, int num) {

int i = 0;

  *n = Calloc1D(llx, sizeof(int));

  for (i = 0; i < num; i++)
    (*n)[xx[i] - 1]++;

}/*FILL_1D_TABLE*/

/* initialize a two-dimensional contingency table and the marginals. */
void fill_2d_table(int *xx, int *yy, int ***n, int **ni, int **nj, int llx,
    int lly, int num) {

int i = 0, j = 0, k = 0;

  *n = (int **) Calloc2D(llx, lly, sizeof(int));
  *ni = (int *) Calloc1D(llx, sizeof(int));
  *nj = (int *) Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    (*n)[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++) {

    (*ni)[i] += (*n)[i][j];
    (*nj)[j] += (*n)[i][j];

  }/*FOR*/

}/*FILL_2D_TABLE*/

/* initialize a three-dimensional contingency table and the marginals. */
void fill_3d_table(int *xx, int *yy, int *zz, int ****n, int ***ni, int ***nj,
    int **nk, int llx, int lly, int llz, int num) {

int i = 0, j = 0, k = 0;

  *n = (int ***) Calloc3D(llz, llx, lly, sizeof(int));
  *ni = (int **) Calloc2D(llz, llx, sizeof(int));
  *nj = (int **) Calloc2D(llz, lly, sizeof(int));
  *nk = (int *) Calloc1D(llz, sizeof(int));

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    (*n)[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      for (k = 0; k < llz; k++) {

        (*ni)[k][i] += (*n)[k][i][j];
        (*nj)[k][j] += (*n)[k][i][j];
        (*nk)[k] += (*n)[k][i][j];

      }/*FOR*/

}/*FILL_3D_TABLE*/

