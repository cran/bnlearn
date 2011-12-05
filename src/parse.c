#include "common.h"

/* find the matching closed brace. */
SEXP match_brace(SEXP lines, SEXP start) {

int depth = 0, open = 0, line_id = INT(start) - 1;
const char *current = NULL;
SEXP stop;

  do {

    /* dereference the current line. */
    current = CHAR(STRING_ELT(lines, line_id));

    /* increment the depth counter if an open brace is found. */
    if (strstr(current, "{")) {

      /* be sure no to exit from the do-while loop until an open curly brace
       * has been spotted. */
      open = 1;
      depth++;

    }/*THEN*/
    /* decrement the depth counter if a closed brace is found. */
    if (strstr(current, "}"))
      depth--;

    /* increment the line id. */
    line_id++;

  } while ((depth > 0) || (open == 0));

  /* allocate and assing the return value. */
  PROTECT(stop = allocVector(INTSXP, 1));
  INT(stop) = line_id;
  UNPROTECT(1);

return stop;

}/*MATCH_BRACE*/

