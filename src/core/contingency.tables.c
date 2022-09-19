#include "../include/rcore.h"
#include "allocations.h"
#include "contingency.tables.h"

/* --------------------- one-dimensional tables ------------------------- */

/* create a one-dimensional contingency table. */
counts1d new_1d_table(int llx) {

counts1d table = { 0 };

  table.llx = llx;
  table.n = Calloc1D(llx, sizeof(int));

  return table;

}/*NEW_1D_TABLE*/

/* initialize a one-dimensional contingency table. */
void fill_1d_table(int *xx, counts1d *table, int num) {

int i = 0, ncomplete = 0;

  /* first fill the counts into the table... */
  for (i = 0; i < num; i++)
    if (xx[i] != NA_INTEGER)
      (*table).n[xx[i] - 1]++;

  /* ... then add them up to count the number of complete observations. */
  for (i = 0; i < (*table).llx; i++)
    ncomplete += (*table).n[i];
  (*table).nobs = ncomplete;

}/*FILL_1D_TABLE*/

/* print a one-dimensional contingency table. */
void print_1d_table(counts1d table) {

  Rprintf("1-dimensional contingency table (%d cells)\n", table.llx);

  for (int i = 0; i < table.llx; i++)
    Rprintf("%d ", table.n[i]);
  Rprintf("\n");

}/*PRINT_1D_TABLE*/

/* free a one-dimensional contingency table. */
void Free1DTAB(counts1d table) {

  Free1D(table.n);

}/*FREE1DTAB*/

/* --------------------- two-dimensional tables ------------------------- */

/* create a two-dimensional contingency table. */
counts2d new_2d_table(int llx, int lly, bool margins) {

counts2d table = { 0 };

  table.llx = llx;
  table.lly = lly;
  table.n = (int **)Calloc2D(llx, lly, sizeof(int));

  if (margins) {

    table.ni = Calloc1D(llx, sizeof(int));
    table.nj = Calloc1D(lly, sizeof(int));

  }/*THEN*/

  return table;

}/*NEW_2D_TABLE*/

/* initialize a two-dimensional contingency table including the marginals. */
void fill_2d_table(int *xx, int *yy, counts2d *table, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    if ((xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*table).n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals if they have been allocated memory for. */
  if ((*table).ni && (*table).nj) {

    for (i = 0; i < (*table).llx; i++)
      for (j = 0; j < (*table).lly; j++) {

      (*table).ni[i] += (*table).n[i][j];
      (*table).nj[j] += (*table).n[i][j];

    }/*FOR*/

    /* compute the number of complete observations. */
    for (i = 0; i < (*table).llx; i++)
      ncomplete += (*table).ni[i];
    (*table).nobs = ncomplete;

  }/*THEN*/
  else {

    /* compute the number of complete observations. */
    for (i = 0; i < (*table).llx; i++)
      for (j = 0; j < (*table).lly; j++)
        ncomplete += (*table).n[i][j];
    (*table).nobs = ncomplete;

  }/*ELSE*/


}/*FILL_2D_TABLE*/

/* same as the above, but zeroes the joint and marginal counts first. */
void refill_2d_table(int *xx, int *yy, counts2d *table, int num) {

  for (int i = 0; i < (*table).llx; i++)
    memset((*table).n[i], '\0', (*table).lly * sizeof(int));
  if ((*table).ni)
    memset((*table).ni, '\0', (*table).llx * sizeof(int));
  if ((*table).nj)
    memset((*table).nj, '\0', (*table).lly * sizeof(int));

  fill_2d_table(xx, yy, table, num);

}/*REFILL_2D_TABLE*/

/* change the dimensions of a two-dimensional contingency table. */
void resize_2d_table(int llx, int lly, counts2d *table) {

  (*table).llx = llx;
  (*table).lly = lly;

}/*RESIZE_2D_TABLE*/

/* print a two-dimensional contingency table. */
void print_2d_table(counts2d table) {

  Rprintf("2-dimensional contingency table (%d x %d cells)\n",
    table.llx, table.lly);
  for (int i = 0; i < table.llx; i++) {

    for (int j = 0; j < table.lly; j++)
      Rprintf("%d ", table.n[i][j]);
    Rprintf("\n");

  }/*FOR*/

}/*PRINT_2D_TABLE*/

/* free a two-dimensional contingency table. */
void Free2DTAB(counts2d table) {

    Free2D(table.n, table.llx);
    Free1D(table.ni);
    Free1D(table.nj);

}/*FREE2DTAB*/

/* --------------------- three-dimensional tables ------------------------ */

/* create a three-dimensional contingency table. */
counts3d new_3d_table(int llx, int lly, int llz) {

counts3d table = { 0 };

  table.llx = llx;
  table.lly = lly;
  table.llz = llz;
  table.n = (int ***) Calloc3D(llz, llx, lly, sizeof(int));
  table.ni = (int **) Calloc2D(llz, llx, sizeof(int));
  table.nj = (int **) Calloc2D(llz, lly, sizeof(int));
  table.nk = (int *) Calloc1D(llz, sizeof(int));

  return table;

}/*NEW_3D_TABLE*/

/* initialize a three-dimensional contingency table and the marginals. */
void fill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num) {

int i = 0, j = 0, k = 0, ncomplete = 0;

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < num; k++)
    if ((zz[k] != NA_INTEGER) && (xx[k] != NA_INTEGER) && (yy[k] != NA_INTEGER))
      (*table).n[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < (*table).llx; i++)
    for (j = 0; j < (*table).lly; j++)
      for (k = 0; k < (*table).llz; k++) {

        (*table).ni[k][i] += (*table).n[k][i][j];
        (*table).nj[k][j] += (*table).n[k][i][j];
        (*table).nk[k] += (*table).n[k][i][j];

      }/*FOR*/

  /* compute the number of complete observations. */
  for (k = 0; k < (*table).llz; k++)
    ncomplete += (*table).nk[k];
  (*table).nobs = ncomplete;

}/*FILL_3D_TABLE*/

/* same as the above, but zeroes the joint and marginal counts first. */
void refill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num) {

  for (int k = 0; k < (*table).llz; k++) {

    for (int i = 0; i < (*table).llx; i++)
      memset((*table).n[k][i], '\0', (*table).lly * sizeof(int));

    memset((*table).ni[k], '\0', (*table).llx * sizeof(int));
    memset((*table).nj[k], '\0', (*table).lly * sizeof(int));

  }/*FOR*/

  memset((*table).nk, '\0', (*table).llz * sizeof(int));

  fill_3d_table(xx, yy, zz, table, num);

}/*REFILL_3D_TABLE*/

/* change the dimensions of a three-dimensional contingency table. */
void resize_3d_table(int llx, int lly, int llz, counts3d *table) {

  (*table).llx = llx;
  (*table).lly = lly;
  (*table).llz = llz;

}/*RESIZE_3D_TABLE*/

/* print a three-dimensional contingency table. */
void print_3d_table(counts3d table) {

  Rprintf("3-dimensional contingency table (%d x %d x %d cells)\n",
    table.llx, table.lly, table.llz);

  for (int k = 0; k < table.llz; k++) {

    Rprintf("[slice %d]", k);

    for (int i = 0; i < table.llx; i++) {

      for (int j = 0; j < table.lly; j++)
        Rprintf("%d ", table.n[k][i][j]);
      Rprintf("\n");

    }/*FOR*/

  }/*FOR*/

}/*PRINT_3D_TABLE*/

/* free a three-dimensional contingency table. */
void Free3DTAB(counts3d table) {

  Free3D(table.n, table.llz, table.llx);
  Free2D(table.ni, table.llz);
  Free2D(table.nj, table.llz);
  Free1D(table.nk);

}/*FREE3DTAB*/

/* --------------------- uniform random table sampling ------------------- */

/* Modified version of the rcont2() function from R. */
static void c_rcont2(int nrow, int ncol, int *nrowt, int *ncolt, int ntotal,
    double *fact, int *jwork, int **matrix) {

int j = 0, l = 0, m = 0, nll = 0, nlm = 0, lsm = 0, lsp = 0;
int ia = 0, ib = 0, ic = 0 , jc = ntotal, id = 0, ie = 0, ii = 0;
double x = 0, y = 0, dummy = 0, sumprb = 0;

  /* Construct random matrix */
  for (j = 0; j < ncol - 1; ++j)
    jwork[j] = ncolt[j];

  for (l = 0; l < nrow - 1; ++l) { /* -----  matrix[ l, * ] ----- */

    ia = nrowt[l];
    ic = jc;
    jc -= ia; /* = n_tot - sum(nr[0:l]) */

    for (m = 0; m < ncol - 1; ++m) {

      id = jwork[m];
      ie = ic;
      ic -= id;
      ib = ie - ia;
      ii = ib - id;

      if (ie == 0) { /* Row [l,] is full, fill rest with zero entries */

        for (j = m; j < ncol - 1; ++j)
          matrix[l][j] = 0;
        ia = 0;
        break;

      }/*FOR*/

      /* Generate pseudo-random number */
      dummy = unif_rand();

      do { /* Outer Loop */

        /* Compute conditional expected value of MATRIX(L, M) */
        nlm = (int)(ia * (id / (double) ie) + 0.5);
        x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id] - fact[ie] - fact[nlm]
                         - fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);
        if (x >= dummy)
          break;
        if (x == 0.) /* MM: I haven't seen this anymore */
          error("rcont2 [%d, %d]: exp underflow to 0; algorithm failure", l, m);

        sumprb = x;
        y = x;
        nll = nlm;

        do {

          /* Increment entry in row L, column M */
          j = (int)((id - nlm) * (double)(ia - nlm));
          lsp = (j == 0);

          if (!lsp) {

            ++nlm;
            x = x * j / ((double) nlm * (ii + nlm));
            sumprb += x;
            if (sumprb >= dummy)
              goto L160;

          }/*THEN*/

          do {

            /* Decrement entry in row L, column M */
            j = (int)(nll * (double)(ii + nll));
            lsm = (j == 0);

            if (!lsm) {

              --nll;
              y = y * j / ((double) (id - nll) * (ia - nll));
              sumprb += y;

              if (sumprb >= dummy) {

                nlm = nll;
                goto L160;

              }/*THEN*/

              if (!lsp)
                break; /* to while (!lsp) */

            }/*THEN*/

          } while (!lsm);

        } while (!lsp);

        dummy = sumprb * unif_rand();

      } while (1);

L160:
      matrix[l][m] = nlm;
      ia -= nlm;
      jwork[m] -= nlm;

    }/*FOR*/

    matrix[l][ncol - 1] = ia;

  }/*FOR*/

  /* Compute entries in last row of MATRIX */
  for (m = 0; m < ncol - 1; ++m)
    matrix[nrow - 1][m] = jwork[m];

  matrix[nrow - 1][ncol - 1] = ib - matrix[nrow - 1][ncol - 2];

}/*C_RCONT2*/

/* generate a random two-dimensional contingency table. */
void rcounts2d(counts2d table, double *fact, int *workspace) {

  c_rcont2(table.llx, table.lly, table.ni, table.nj, table.nobs, fact,
        workspace, table.n);

}/*RCOUNTS2D*/

/* generate a random three-dimensional contingency table. */
void rcounts3d(counts3d table, double *fact, int *workspace) {

  for (int k = 0; k < table.llz; k++)
    c_rcont2(table.llx, table.lly, table.ni[k], table.nj[k], table.nk[k],
      fact, workspace, table.n[k]);

}/*RCOUNTS3D*/

