#include "common.h"

/* macro for the number of levels of the j-th node. */
#define BNLEARN_NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(data, j), R_LevelsSymbol))

/* get the number of parameters of a single node (discrete case). */
SEXP nparams_dnode(SEXP graph, SEXP node, SEXP data, SEXP real) {

int i = 0, j = 0, nlevels = 1;
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

      nlevels *= BNLEARN_NLEVELS(i) - 1 * INT(real);

    }/*THEN*/

  }/*FOR*/

  PROTECT(temp = allocVector(INTSXP, 1));
  INT(temp) = nlevels;
  UNPROTECT(1);

return temp;

}/*NPARAMS_DNODE*/

/* get the number of parameters of a single node (discrete case). */
SEXP nparams_gnode(SEXP graph, SEXP node) {

char *name = (char *)CHAR(STRING_ELT(node, 0));
SEXP temp, result;

  temp = getListElement(graph, "nodes");
  temp = getListElement(temp, name);
  temp = getListElement(temp, "parents");

  PROTECT(result = allocVector(INTSXP, 1));
  INT(result) = LENGTH(temp) + 1;
  UNPROTECT(1);

  return result;

}/*NPARAMS_GNODE*/

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

/* convert a set of neighbourhoods into an arc set. */
SEXP nbr2arcs(SEXP nbr) {

int i = 0, j = 0, k = 0, narcs = 0;
int length_names = 0;
SEXP arcs, dimnames, colnames, temp, names;

  /* get the names of the nodes. */
  names = getAttrib(nbr, R_NamesSymbol);
  length_names = LENGTH(names);

  /* scan the structure to determine the number of arcs.  */
  for (i = 0; i < length_names; i++) {

    /* get the entry for the neighbours of the node.*/
    temp = getListElement(nbr, (char *)CHAR(STRING_ELT(names, i)));
    temp = getListElement(temp, "nbr");

    narcs += LENGTH(temp);

  }/*FOR*/

  /* allocate colnames. */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("from"));
  SET_STRING_ELT(colnames, 1, mkChar("to"));
  SET_VECTOR_ELT(dimnames, 1, colnames);

  /* if there are no arcs, return an empty arc set. */
  if (narcs == 0) {

    /* allocate an empty arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, 0, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

    UNPROTECT(3);

    return arcs;

  }/*THEN*/
  else {

    /* allocate the arc set. */
    PROTECT(arcs = allocMatrix(STRSXP, narcs, 2));
    /* set the column names. */
    setAttrib(arcs, R_DimNamesSymbol, dimnames);

  }/*ELSE*/

  /* rescan the structure to build the arc set. */
  for (i = 0; i < length_names; i++) {

    /* get the entry for the neighbours of the node.*/
    temp = getListElement(nbr, (char *)CHAR(STRING_ELT(names, i)));
    temp = getListElement(temp, "nbr");

    for (j = 0; j < LENGTH(temp); j++) {

      SET_STRING_ELT(arcs, k, STRING_ELT(names, i));
      SET_STRING_ELT(arcs, k + 1 * narcs , STRING_ELT(temp, j));
      k++;

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(3);

return arcs;

}/*NBR2ARCS*/

