
#include "common.h"

/* unconditional mutual information, to be used in C code. */
double c_mi(int *xx, int *llx, int *yy, int *lly, int *num) {

int i = 0, j = 0, k = 0;
int  **n = NULL, *ni = NULL, *nj = NULL;
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(*llx, *lly);
  ni = alloc1dcont(*llx);
  nj = alloc1dcont(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) 
      res += MI_PART(n[i][j], ni[i], nj[j], *num);

  return (res)/(*num);

}/*C_MI*/

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length) {

int *llx = INTEGER(lx), *lly = INTEGER(ly), *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y);
SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));
  NUM(result) = c_mi(xx, llx, yy, lly, num);
  UNPROTECT(1);

  return result;

}/*MI*/

/* conditional mutual information, to be used in C code. */
double c_cmi(int *xx, int *llx, int *yy, int *lly, int *zz, int *llz, int *num) {

int i = 0, j = 0, k = 0; 
int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dcont(*llx, *lly, *llz);
  ni = alloc2dcont(*llx, *llz);
  nj = alloc2dcont(*lly, *llz);
  nk = alloc1dcont(*llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1][zz[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) {

        ni[i][k] += n[i][j][k];
        nj[j][k] += n[i][j][k];
        nk[k] += n[i][j][k];

      }/*FOR*/

  /* compute the conditional mutual information from the joint and
     marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) 
        res += MI_PART(n[i][j][k], ni[i][k], nj[j][k], nk[k]);

  res = res/(*num);

  return res;

}/*C_CMI*/

/* conditional mutual information, to be used for the asymptotic test. */
SEXP cmi(SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length) {

int *llx = INTEGER(lx), *lly = INTEGER(ly), *llz = INTEGER(lz);
int *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
SEXP result;

  /* allocate and initialize result to zero. */
  PROTECT(result = allocVector(REALSXP, 1));
  NUM(result) = c_cmi(xx, llx, yy, lly, zz, llz, num);
  UNPROTECT(1);

  return result;

}/*CMI*/

/* unconditional Gaussian mutual information, to be used in C code. */
double c_mig(double *xx, double *yy, int *num) {

double cor = c_fast_cor(xx, yy, num);

  return - 0.5 * log(1 - cor * cor);

}/*C_MIG*/

/* unconditional Gaussian mutual information, to be used in the asymptotic test. */
SEXP mig(SEXP x, SEXP y, SEXP length) {

double *xx = REAL(x), *yy = REAL(y);
int *num = INTEGER(length);
SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));
  NUM(result) = c_mig(xx, yy, num);
  UNPROTECT(1);

  return result;

}/*MIG*/
