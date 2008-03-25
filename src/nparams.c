#include <R.h>
#include <Rinternals.h>

#define BNLEARN_NAMES(x) \
  LENGTH(getAttrib(x, R_NamesSymbol))

#define BNLEARN_NAME(x, i) \
  CHAR(STRING_ELT(getAttrib(x, R_NamesSymbol), i))

/* macro for the number of levels of the j-th node. */
#define BNLEARN_NLEVELS(j) \
  LENGTH(getAttrib(VECTOR_ELT(data, j), R_LevelsSymbol))

SEXP nparams(SEXP graph, SEXP node, SEXP data, SEXP real) {

  int i = 0, j = 0;
  int nlevels = 1;

  SEXP temp;

  /* loop over the list names. */
  for (i = 0; i < BNLEARN_NAMES(graph); i++) {

    /* the element I need is called "nodes". */
    if (!strcmp(BNLEARN_NAME(graph, i), "nodes")) {

      /* "nodes" is a list itself, set a SEXP for later use. */
      temp = VECTOR_ELT(graph, i);

      break;

    }/*THEN*/

  }/*FOR*/

  /* loop over the nodes of the graph. */
  for (i = 0; i < BNLEARN_NAMES(temp); i++) {

    /* the element I need is stored in SEXP node. */
    if (!strcmp(BNLEARN_NAME(temp, i), CHAR(STRING_ELT(node, 0)))) {

      /* each node is a list, set a SEXP for later use. */
      temp = VECTOR_ELT(temp, i);

      break;

    }/*THEN*/

  }/*FOR*/

  /* loop over node cached information. */
  for (i = 0; i < BNLEARN_NAMES(temp); i++) {

    /* the element I need is called "parents".. */
    if (!strcmp(BNLEARN_NAME(temp, i), "parents")) {

      /* the parents are stored in a character vector, set a SEXP
         for later use. */
      temp = VECTOR_ELT(temp, i);

      break;

    }/*THEN*/

  }/*FOR*/

  /* sum (multiply, actually) up the levels. */
  for (i = 0; i < BNLEARN_NAMES(data); i++) {

    for (j = 0; j < LENGTH(temp); j++) {

      /* this is a parent. */
      if (!strcmp(BNLEARN_NAME(data, i), CHAR(STRING_ELT(temp, j)))) {

        nlevels *= BNLEARN_NLEVELS(i);

      }/*THEN*/

    }/*FOR*/

    /* this is the node. */
    if (!strcmp(BNLEARN_NAME(data, i), CHAR(STRING_ELT(node, 0)))) {

      nlevels *= BNLEARN_NLEVELS(i) - 1 * INTEGER(real)[0];

    }/*THEN*/

  }/*FOR*/

  PROTECT(temp = allocVector(INTSXP, 1));
  INTEGER(temp)[0] = nlevels;
  UNPROTECT(1);

return temp;

}/*NPARAMS*/

