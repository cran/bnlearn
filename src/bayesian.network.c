#include <R.h>
#include <Rinternals.h>

/* macro for the number of levels of the j-th node. */
#define BNLEARN_NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(data, j), R_LevelsSymbol))

/* get the list element named str, or return NULL. */
SEXP getListElement(SEXP list, char *str) {

  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i = 0;

  for (i = 0; i < length(list); i++) {

    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {

      elmt = VECTOR_ELT(list, i);
      break;

    }/*THEN*/

  }/*FOR*/

return elmt;

}/*GETLISTELEMENT*/

/* get the number of parameters of a single node. */
SEXP nparams(SEXP graph, SEXP node, SEXP data, SEXP real) {

  int i = 0, j = 0;
  int nlevels = 1;
  int length_names = 0, length_nodes = 0;

  SEXP temp, names;

  /* get the entry for the parents of the node.*/
  temp = getListElement(graph, "nodes");
  temp = getListElement(temp, (char *)CHAR(STRING_ELT(node, 0)));
  temp = getListElement(temp, "parents");

  /* get the column names from the data set and the length of the
       relevant vectors. */
  names = getAttrib(data, R_NamesSymbol);
  length_names = LENGTH(names);
  length_nodes = LENGTH(temp);

  /* sum (multiply, actually) up the levels. */
  for (i = 0; i < length_names; i++) {

    for (j = 0; j < length_nodes; j++) {

      /* this is a parent. */
      if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(temp, j)))) {

        nlevels *= BNLEARN_NLEVELS(i);

      }/*THEN*/

    }/*FOR*/

    /* this is the node. */
    if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(node, 0)))) {

      nlevels *= BNLEARN_NLEVELS(i) - 1 * INTEGER(real)[0];

    }/*THEN*/

  }/*FOR*/

  PROTECT(temp = allocVector(INTSXP, 1));
  INTEGER(temp)[0] = nlevels;
  UNPROTECT(1);

return temp;

}/*NPARAMS*/

/* schedule the children of the current nodes in a breadth-first search. */
SEXP schedule_children(SEXP graph, SEXP nodes) {

  SEXP temp, names, children, result;
  int count = 0, i = 0, j = 0, k = 0, l = 0;
  int length_nodes = 0, length_names = 0;

  /* get the nodes' structures and their names. */
  temp = getListElement(graph, "nodes");
  names = getAttrib(temp, R_NamesSymbol);

  /* compute the length of both vectors, once.*/
  length_nodes = LENGTH(nodes);
  length_names = LENGTH(names);

  /* count how many they are. */
  for (i = 0; i < length_names; i++) {

    for (j = 0; j < length_nodes; j ++) {

      if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(nodes, j)))) {

        children = getListElement(temp, (char *)CHAR(STRING_ELT(nodes, j)));
        children = getListElement(children, "children");
        count += LENGTH(children);

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  /* allocate and protect the result. */
  PROTECT(result = allocVector(STRSXP, count));

  if (count == 0) {

    UNPROTECT(1);
    return result;

  }/*THEN*/

  /* fill in the result vector. */
  for (i = 0; i < length_names; i++) {

    for (j = 0; j < length_nodes; j ++) {

      if (!strcmp(CHAR(STRING_ELT(names, i)), CHAR(STRING_ELT(nodes, j)))) {

        children = getListElement(temp, (char *)CHAR(STRING_ELT(nodes, j)));
        children = getListElement(children, "children");

        /* if there are no children, skip to the next node. */
        if (isNull(children)) continue;

        for (k = 0; k < LENGTH(children); k++)
          SET_STRING_ELT(result, l++, STRING_ELT(children, k));

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

 UNPROTECT(1);

 return result;

}/*SCHEDULE_CHILDREN*/

