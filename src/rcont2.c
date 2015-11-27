/* Modified version of the rcont2() function from R. */
#include "include/rcore.h"

void c_rcont2(int nrow, int ncol, int *nrowt, int *ncolt, int ntotal,
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

            R_CheckUserInterrupt();

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

