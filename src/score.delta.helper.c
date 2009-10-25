
#include "common.h"

SEXP score_delta_helper(SEXP net, SEXP arc, SEXP operator) {

int i = 0, k = 0;
char *from = (char *) CHAR(STRING_ELT(arc, 0));
char *to = (char *) CHAR(STRING_ELT(arc, 1));
char *op = (char *) CHAR(STRING_ELT(operator, 0));
SEXP fake, nodes, parents_from, parents_to;
SEXP start, cur, temp, _to, _nodes, _parents;

  /* allocate the return value. */
  PROTECT(fake = allocVector(VECSXP, 1));

  /* allocate and initialize the names of the elements of the fake network. */
  PROTECT(_parents = allocVector(STRSXP, 1));
  SET_STRING_ELT(_parents, 0, mkChar("parents"));
  PROTECT(_nodes = allocVector(STRSXP, 1));
  SET_STRING_ELT(_nodes, 0, mkChar("nodes"));

  start = getListElement(net, "nodes");

  if (!strcmp(op, "set")) {

    /* adding the arc to the graph. */
    PROTECT(nodes = allocVector(VECSXP, 1));
    PROTECT(temp = allocVector(VECSXP, 1));
    PROTECT(_to = allocVector(STRSXP, 1));
    SET_STRING_ELT(_to, 0, mkChar(to));

    setAttrib(nodes, R_NamesSymbol, _to);
    setAttrib(temp, R_NamesSymbol, _parents);

    cur = getListElement(start, to);
    cur = getListElement(cur, "parents");

    PROTECT(parents_to = allocVector(STRSXP, LENGTH(cur) + 1));

    /* the new parents are the old ones plus the other node incident on
     * the arc. */
    for (i = 0; i < LENGTH(cur); i++)
      SET_STRING_ELT(parents_to, i, STRING_ELT(cur, i));
    SET_STRING_ELT(parents_to, LENGTH(cur), STRING_ELT(arc, 0));

    SET_VECTOR_ELT(temp, 0, parents_to);
    SET_VECTOR_ELT(nodes, 0, temp);

    UNPROTECT(4);

  }/*THEN*/
  else if (!strcmp(op, "drop")) {

    /* dropping the arc from the graph. */
    PROTECT(nodes = allocVector(VECSXP, 1));
    PROTECT(temp = allocVector(VECSXP, 1));
    PROTECT(_to = allocVector(STRSXP, 1));
    SET_STRING_ELT(_to, 0, mkChar(to));
    setAttrib(nodes, R_NamesSymbol, _to);
    setAttrib(temp, R_NamesSymbol, _parents);

    cur = getListElement(start, to);
    cur = getListElement(cur, "parents");

    PROTECT(parents_to = allocVector(STRSXP, LENGTH(cur) - 1));

    /* the new parents are the old ones except the other node incident on
     * the arc. */
    for (i = 0, k = 0; i < LENGTH(cur); i++)
      if (strcmp(CHAR(STRING_ELT(cur, i)), from) != 0)
        SET_STRING_ELT(parents_to, k++, STRING_ELT(cur, i));

    SET_VECTOR_ELT(temp, 0, parents_to);
    SET_VECTOR_ELT(nodes, 0, temp);

    UNPROTECT(4);

  }/*THEN*/
  else {

    /* reversing the arc in the graph. */
    PROTECT(nodes = allocVector(VECSXP, 2));
    PROTECT(temp = allocVector(VECSXP, 1));
    setAttrib(nodes, R_NamesSymbol, arc);
    setAttrib(temp, R_NamesSymbol, _parents);

    /* add "to" to the parents of "from". */
    cur = getListElement(start, from);
    cur = getListElement(cur, "parents");

    PROTECT(parents_from = allocVector(STRSXP, LENGTH(cur) + 1));

    for (i = 0; i < LENGTH(cur); i++)
      SET_STRING_ELT(parents_from, i, STRING_ELT(cur, i));
    SET_STRING_ELT(parents_from, LENGTH(cur), STRING_ELT(arc, 1));

    SET_VECTOR_ELT(temp, 0, parents_from);
    SET_VECTOR_ELT(nodes, 0, duplicate(temp));

    /* remove "from" from the parents of "to". */
    cur = getListElement(start, to);
    cur = getListElement(cur, "parents");

    PROTECT(parents_to = allocVector(STRSXP, LENGTH(cur) - 1));

    for (i = 0, k = 0; i < LENGTH(cur); i++)
      if (strcmp(CHAR(STRING_ELT(cur, i)), from) != 0)
        SET_STRING_ELT(parents_to, k++, STRING_ELT(cur, i));

    SET_VECTOR_ELT(temp, 0, parents_to);
    SET_VECTOR_ELT(nodes, 1, temp);

    UNPROTECT(4);

  }/*ELSE*/

  /* save the fabricated nodes' strctures in the return value. */
  SET_VECTOR_ELT(fake, 0, nodes);
  setAttrib(fake, R_NamesSymbol, _nodes);

  UNPROTECT(3);

  return fake;

}/*SCORE_DELTA_HELPER*/

