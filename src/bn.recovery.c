#include "include/rcore.h"
#include "include/matrix.h"

/* check neighbourhood sets and markov blanets for consistency.. */
SEXP bn_recovery(SEXP bn, SEXP mb, SEXP filter, SEXP debug) {

int i = 0, j = 0, k = 0, n = 0, counter = 0, *flt = INTEGER(filter);
short int *checklist = NULL, err = 0;
bool debugging = isTRUE(debug), checkmb = isTRUE(mb);
SEXP temp, temp2, nodes, elnames = NULL, fixed;

  /* get the names of the nodes. */
  PROTECT(nodes = getAttrib(bn, R_NamesSymbol));
  n = length(nodes);

  /* allocate and initialize the checklist. */
  checklist = Calloc1D(UPTRI_MATRIX(n), sizeof(short int));

  if (debugging) {

    Rprintf("----------------------------------------------------------------\n");

    if (checkmb)
      Rprintf("* checking consistency of markov blankets.\n");
    else
      Rprintf("* checking consistency of neighbourhood sets.\n");

   }/*THEN*/

  /* scan the structure to determine the number of arcs.  */
  for (i = 0; i < n; i++) {

     if (debugging)
       Rprintf("  > checking node %s.\n",  NODE(i));

    /* get the entry for the (neighbours|elements of the markov blanket)
       of the node.*/
    temp = getListElement(bn, (char *)NODE(i));
    if (!checkmb)
      temp = getListElement(temp, "nbr");

    /* check each element of the array and identify which variable it
       corresponds to. */
    for (j = 0; j < length(temp); j++) {

      for (k = 0; k < n; k++) {

        /* increment the right element of checklist. */
        if (!strcmp(NODE(k), (char *)CHAR(STRING_ELT(temp, j))))
          checklist[UPTRI(i + 1, k + 1, n)]++;

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

      if ((checklist[UPTRI(i + 1, j + 1, n)] != 0) &&
          (checklist[UPTRI(i + 1, j + 1, n)] != 2)) {

        if (debugging) {

          if (checkmb)
            Rprintf("@ asymmetry in the markov blankets for %s and %s.\n",
              NODE(i), NODE(j));
          else
            Rprintf("@ asymmetry in the neighbourhood sets for %s and %s.\n",
              NODE(i), NODE(j));

        }/*THEN*/

        err = 1;

      }/*THEN*/

    }/*FOR*/

  /* no need to go on if the (neighbourhoods|markov blankets) are symmetric. */
  if (!err) {

    Free1D(checklist);
    UNPROTECT(1);

    return bn;

  }/*THEN*/

  /* build a correct structure to return. */
  PROTECT(fixed = allocVector(VECSXP, n));
  setAttrib(fixed, R_NamesSymbol, nodes);

  /* allocate colnames. */
  if (!checkmb)
    PROTECT(elnames = mkStringVec(2, "mb", "nbr"));

  for (i = 0; i < n; i++) {

    if (!checkmb) {

      /* allocate the "mb" and "nbr" elements of the node. */
      PROTECT(temp = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(fixed, i, temp);
      setAttrib(temp, R_NamesSymbol, elnames);

      /* copy the "mb" part from the old structure. */
      temp2 = getListElement(bn, (char *)NODE(i));
      temp2 = getListElement(temp2, "mb");
      SET_VECTOR_ELT(temp, 0, temp2);

    }/*THEN*/

    /* rescan the checklist. */
    for (j = 0; j < n; j++)
      if (checklist[UPTRI(i + 1, j + 1, n)] >= *flt)
        if (i != j)
          counter++;

    /* allocate and fill the "nbr" element. */
    PROTECT(temp2 = allocVector(STRSXP, counter));

    for (j = 0; j < n; j++)
      if (checklist[UPTRI(i + 1, j + 1, n)] == *flt)
        if (i != j)
          SET_STRING_ELT(temp2, --counter, STRING_ELT(nodes, j));

    if (checkmb) {

      SET_VECTOR_ELT(fixed, i, temp2);
      UNPROTECT(1);

    }/*THEN*/
    else {

      SET_VECTOR_ELT(temp, 1, temp2);
      UNPROTECT(2);

    }/*ELSE*/

  }/*FOR*/

  if (checkmb)
    UNPROTECT(2);
  else
    UNPROTECT(3);

  Free1D(checklist);

  return fixed;

}/*BN_RECOVERY*/

