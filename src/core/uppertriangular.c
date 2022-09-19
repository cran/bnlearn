#include "../include/rcore.h"
#include "allocations.h"
#include "uppertriangular.h"

uppertriangular new_uppertriangular(int dim) {

uppertriangular sym = { 0 };

  sym.dim = dim;
  sym.mat = Calloc1D(dim * (dim - 1) / 2, sizeof(double));

  return sym;

}/*NEW_UPPERTRIANGULAR*/

void FreeUPPERTRIANGULAR(uppertriangular sym) {

  Free1D(sym.names);
  Free1D(sym.mat);

}/*FREEUPPERTRIANGULAR*/

int uppertriangular_size(uppertriangular sym) {

  return sym.dim * (sym.dim - 1) / 2;

}/*UPPTERTRIANGULAR_SIZE*/

void uppertriangular_copy_names(uppertriangular *sym, const char **names) {

int i = 0;

  (*sym).names = Calloc1D((*sym).dim, sizeof(char *));

  for (i = 0; i < (*sym).dim; i++)
    (*sym).names[i] = names[i];

}/*UPPERTRIANGULAR_COPY_NAMES*/

/* inverse function of the UPTRI3() macro. */
void INV_UPTRI3(int x, int n, int *res) {

int c = 0, r = 0, cn = n - 1;

  for (r = 0; r < n; r++) {

    if (x < cn) {

      c = n - (cn - x);
      break;

    }/*THEN*/
    else {

      cn += n - (r + 2);

    }/*ELSE*/

  }/*FOR*/

  res[0] = r;
  res[1] = c;

}/*INV_UPTRI3*/

