#include "include/rcore.h"
#include "include/allocations.h"

/* initialize a two-dimensional contingency table and the marginals. */
void fill_2d_table(int *xx, int *yy, int ***n, int **ni, int **nj, int llx,
    int lly, int num) {

int i = 0, j = 0, k = 0;

  *n = alloc2dcont(llx, lly);
  *ni = alloc1dcont(llx);
  *nj = alloc1dcont(lly);

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

  *n = alloc3dcont(llz, llx, lly);
  *ni = alloc2dcont(llz, llx);
  *nj = alloc2dcont(llz, lly);
  *nk = alloc1dcont(llz);

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

