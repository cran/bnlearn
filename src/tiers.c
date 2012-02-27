#include "common.h"

SEXP tiers(SEXP nodes, SEXP debug) {

int i = 0, j = 0, k = 0, narcs = 0, nnodes = 0, ntiers = LENGTH(nodes);
int *tier_size = NULL, *debuglevel = LOGICAL(debug), tier_start = 0, cur = 0;
SEXP flattened, blacklist, temp;

  /* allocate the counters for tiers' sizes.*/
  tier_size = alloc1dcont(ntiers);

  if (!isString(nodes)) {

    /* "node" is a list, each tier is an element. */
    for (i = ntiers - 1; i >= 0; i--) {

      temp = VECTOR_ELT(nodes, i);
      tier_size[i] = LENGTH(temp);  
      nnodes += tier_size[i];
      narcs += (nnodes - tier_size[i]) * tier_size[i];

    }/*FOR*/

    /* flatten the tiers to keep manipulation later on simple. */
    PROTECT(flattened = allocVector(STRSXP, nnodes));

    for (i = 0, k = 0; i < ntiers; i++) {

      temp = VECTOR_ELT(nodes, i);

      for (j = 0; j < tier_size[i]; j++)
        SET_STRING_ELT(flattened, k++, STRING_ELT(temp, j));

    }/*FOR*/

  }/*THEN*/
  else {

    /* "node" is a character vector, which means that each node is in its own tier
     * and that there is no need to flatted it. */
    flattened = nodes; 
    nnodes = LENGTH(nodes);
    for (i = 0; i < ntiers; i++)
      tier_size[i] = 1;

    /* the blacklist is the one resulting from a complete node ordering. */
    narcs = ntiers * (ntiers - 1) / 2;

  }/*ELSE*/

  /* allocate the return value. */
  PROTECT(blacklist = allocMatrix(STRSXP, narcs, 2));

  for (k = 0, i = 0; k < nnodes; k++) {

    temp = STRING_ELT(flattened, k);

    if (*debuglevel > 0)
      Rprintf("* current node is %s in tier %d.\n", CHAR(temp), i + 1);

    for (j = tier_start + tier_size[i]; j < nnodes; j++) {

      if (*debuglevel)
        Rprintf("  > blacklisting %s -> %s\n", CHAR(STRING_ELT(flattened, j)), CHAR(temp));

      SET_STRING_ELT(blacklist, cur, STRING_ELT(flattened, j));
      SET_STRING_ELT(blacklist, cur + narcs, temp);
      cur++;

    }/*FOR*/

    if (k >= tier_start + tier_size[i] - 1)
      tier_start += tier_size[i++];

    if (i == ntiers)
      break;

  }

  /* allocate, initialize and set the column names. */
  finalize_arcs(blacklist);

  if (!isString(nodes))
    UNPROTECT(2);
  else 
    UNPROTECT(1);

  return blacklist;

}/*TIERS*/

