#include "include/rcore.h"
#include "include/graph.h"
#include "include/globals.h"

/* sanitize an arc set using the conditional Gaussian assumptions. */
SEXP arcs_cg_assumptions(SEXP arcs, SEXP nodes, SEXP data) {

int i = 0, k = 0, narcs = length(arcs)/2, nnodes = length(data), ndrop = 0;
int *t = NULL, *und = NULL, *vartype = NULL, *arc_ok = NULL;
SEXP result, try, undirected;

  /* preallocate an array for the typeiable types, to be used to avoid
   * duplicate VECTOR_ELT() calls and thus probe each node just once. */
  vartype = Calloc1D(nnodes, sizeof(int));
  arc_ok = Calloc1D(narcs, sizeof(int));
  /* match the nodes in the arcs to the data. */
  PROTECT(try = match(nodes, arcs, 0));
  t = INTEGER(try);

  /* cache the variables' types. */
  for (i = 0; i < narcs; i++) {

    if (vartype[t[i] - 1] == 0)
      vartype[t[i] - 1] = TYPEOF(VECTOR_ELT(data, t[i] - 1));
    if (vartype[t[i + narcs] - 1] == 0)
      vartype[t[i + narcs] - 1] = TYPEOF(VECTOR_ELT(data, t[i + narcs] - 1));

  }/*FOR*/

  /* check which arcs are undirected, they are sanitized differently. */
  PROTECT(undirected = which_undirected(arcs, nodes));
  und = INTEGER(undirected);

  /* check variables types in each arc are ok. */
  for (i = 0; i < narcs; i++) {

    arc_ok[i] = !((vartype[t[i] - 1] == REALSXP) &&
                  (vartype[t[i + narcs] - 1] == INTSXP));

    /* nothing to do here, move to the next arc. */
    if (arc_ok[i])
      continue;

    /* if the arc violates the assumptions and it's directed, throw an error. */
    if (!und[i]) {

      Free1D(vartype);
      Free1D(arc_ok);
      UNPROTECT(2);

      error("arc %s -> %s violates the assumptions of the model.",
        NODE(t[i] - 1), NODE(t[i + narcs] - 1));

    }/*THEN*/

    /* if the arc is undirected, remove the offending direction. */
    warning("the direction %s -> %s of %s - %s violates the assumptions of the model and will be ignored.",
      NODE(t[i] - 1), NODE(t[i + narcs] - 1), NODE(t[i] - 1), NODE(t[i + narcs] - 1));

    arc_ok[i] = FALSE;
    ndrop++;

  }/*FOR*/

  UNPROTECT(2);

  PROTECT(result = allocMatrix(STRSXP, narcs - ndrop, 2));

  for (i = 0; i < narcs; i++) {

    if (!arc_ok[i])
      continue;

    SET_STRING_ELT(result, k, STRING_ELT(arcs, i));
    SET_STRING_ELT(result, k + narcs - ndrop, STRING_ELT(arcs, i + narcs));
    k++;

  }/*FOR*/

  UNPROTECT(1);
  Free1D(vartype);
  Free1D(arc_ok);

  return arcs;

}/*ARCS_CG_ASSUMPTIONS*/

/* enumerate arcs that violate conditional Gaussian assumptions. */
SEXP cg_banned_arcs(SEXP nodes, SEXP variables) {

int i = 0, j = 0, k = 0, ndp = 0, nnodes = length(nodes), *type = NULL;
const char *class = CHAR(STRING_ELT(getAttrib(variables, R_ClassSymbol), 0));
const char *var_class = NULL;
SEXP split, dpar, gpar, arcs;

  /* cache the variables' types from data frames and fitted networks. */
  type = Calloc1D(nnodes, sizeof(int));

  if (strcmp(class, "data.frame") == 0) {

    for (i = 0; i < nnodes; i++) {

      type[i] = TYPEOF(VECTOR_ELT(variables, i));
      ndp += (type[i] == INTSXP);

    }/*FOR*/

  }/*THEN*/
  else {

    for (i = 0; i < nnodes; i++) {

      var_class = CHAR(STRING_ELT(getAttrib(VECTOR_ELT(variables, i), R_ClassSymbol), 0));

      if (strcmp(var_class, "bn.fit.dnode") == 0) {

        type[i] = INTSXP;
        ndp++;

      }/*THEN*/
      else {

        type[i] = REALSXP;

      }/*ELSE*/

    }/*FOR*/

  }/*ELSE*/

  PROTECT(split = allocVector(VECSXP, 2));
  PROTECT(dpar = allocVector(STRSXP, ndp));
  PROTECT(gpar = allocVector(STRSXP, nnodes - ndp));
  SET_VECTOR_ELT(split, 0, dpar);
  SET_VECTOR_ELT(split, 1, gpar);

  for (i = 0, j = 0, k = 0; i < nnodes; i++)
    if (type[i] == INTSXP)
      SET_STRING_ELT(dpar, j++, STRING_ELT(nodes, i));
    else
      SET_STRING_ELT(gpar, k++, STRING_ELT(nodes, i));

  arcs = tiers(split, FALSESEXP);

  Free1D(type);
  UNPROTECT(3);

  return arcs;

}/*CG_BANNED_ARCS*/
