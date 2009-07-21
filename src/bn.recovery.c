#include "common.h"

#define NODE(i) (char *)CHAR(STRING_ELT(names, i))

/* check neighbourhood sets and markov blanets for consistency.. */
SEXP bn_recovery(SEXP bn, SEXP strict, SEXP mb, SEXP debug) {

  int i = 0, j = 0, k = 0, n = 0;
  short int *checklist;
  int counter = 0;
  short int err = 0;

  SEXP temp, temp2, names, elnames = NULL, fixed;

  /* get the names of the nodes. */
  names = getAttrib(bn, R_NamesSymbol);
  n = LENGTH(names);

  /* allocate and initialize the checklist. */
  checklist = allocstatus(UPTRI(n, n, n));

  if (isTRUE(debug)) {

    Rprintf("----------------------------------------------------------------\n");

    if (isTRUE(mb))
      Rprintf("* checking consistency of markov blankets.\n");
    else
      Rprintf("* checking consistency of neighbourhood sets.\n");

   }/*THEN*/

  /* scan the structure to determine the number of arcs.  */
  for (i = 0; i < n; i++) {

     if (isTRUE(debug))
       Rprintf("  > checking node %s.\n",  NODE(i));

    /* get the entry for the (neighbours|elements of the markov blanket)
       of the node.*/
    temp = getListElement(bn, NODE(i));
    if (!isTRUE(mb))
      temp = getListElement(temp, "nbr");

    /* check each element of the array and identify which variable it
       corresponds to. */
    for (j = 0; j < LENGTH(temp); j++) {

      for (k = 0; k < n; k++) {

        /* increment the right element of checklist. */
        if (!strcmp(NODE(k), (char *)CHAR(STRING_ELT(temp, j))))
          checklist[UPTRI(i + 1, k + 1, n) - 1]++;

      }/*FOR*/

    }/*FOR*/

  }/*FOR*/

  /* if A is a neighbour of B, B is a neighbour of A; therefore each entry in
   * the checklist array must be equal to either zero (if the corresponding
   * nodes are not neighbours) or two (if the corresponding nodes are neighbours).
   * Any other value (typically one) is caused by an incorrect (i.e. asymmetric)
   * neighbourhood structure. The same logic holds for the markov blankets. */
  for (i = 0; i < n; i++)
    for (j = i; j < n; j++) {

      if ((checklist[UPTRI(i + 1, j + 1, n) - 1] != 0) &&
          (checklist[UPTRI(i + 1, j + 1, n) - 1] != 2)) {

        if (isTRUE(debug)) {

          if (isTRUE(mb))
            Rprintf("@ asymmetry in the markov blankets for %s and %s.\n",
              NODE(i), NODE(j));
          else
            Rprintf("@ asymmetry in the neighbourhood sets for %s and %s.\n",
              NODE(i), NODE(j));

        }/*THEN*/

        err = 1;

      }/*THEN*/

    }/*FOR*/

  /* no need to go on if the (neighbourhood sets|markov blankets) are symmetric;
   * otherwise throw either an error or a warning according to the value of the
   * strict parameter. */
  if (!err) {

    return bn;

  }/*THEN*/
  else if (isTRUE(strict)) {

    if (isTRUE(mb))
      error("markov blankets are not symmetric.\n");
    else
      error("neighbourhood sets are not symmetric.\n");

  }/*THEN*/
  else {

    if (isTRUE(mb))
      warning("markov blankets are not symmetric.\n");
    else
      warning("neighbourhood sets are not symmetric.\n");

  }/*ELSE*/

  /* build a correct structure to return. */
  PROTECT(fixed = allocVector(VECSXP, n));
  setAttrib(fixed, R_NamesSymbol, names);

  if (!isTRUE(mb)) {

    /* allocate colnames. */
    PROTECT(elnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(elnames, 0, mkChar("mb"));
    SET_STRING_ELT(elnames, 1, mkChar("nbr"));

  }/*THEN*/

  for (i = 0; i < n; i++) {

    if (!isTRUE(mb)) {

      /* allocate the "mb" and "nbr" elements of the node. */
      PROTECT(temp = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(fixed, i, temp);
      setAttrib(temp, R_NamesSymbol, elnames);

      /* copy the "mb" part from the old structure. */
      temp2 = getListElement(bn, NODE(i));
      temp2 = getListElement(temp2, "mb");
      SET_VECTOR_ELT(temp, 0, temp2);

    }/*THEN*/

    /* rescan the checklist. */
    for (j = 0; j < n; j++)
      if (checklist[UPTRI(i + 1, j + 1, n) - 1] == 2)
        if (i != j)
          counter++;

    /* allocate and fill the "nbr" element. */
    PROTECT(temp2 = allocVector(STRSXP, counter));

    for (j = 0; j < n; j++)
      if (checklist[UPTRI(i + 1, j + 1, n) - 1] == 2)
        if (i != j)
          SET_STRING_ELT(temp2, --counter, STRING_ELT(names, j));

    if (isTRUE(mb)) {

      SET_VECTOR_ELT(fixed, i, temp2);
      UNPROTECT(1);

    }/*THEN*/
    else {

      SET_VECTOR_ELT(temp, 1, temp2);
      UNPROTECT(2);

    }/*ELSE*/

  }/*FOR*/

  if (isTRUE(mb))
    UNPROTECT(1);
  else
    UNPROTECT(2);

return fixed;

}/*BN_RECOVERY*/

