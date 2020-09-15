#include "include/rcore.h"
#include "include/data.frame.h"
#include "include/sets.h"
#include "include/contingency.tables.h"

static int deff_node_root(int *n, int llx) {

int i = 0, nonzero = 0;

  for (i = 0; i < llx; i++)
    nonzero += (n[i] > 0);

  return nonzero;

}/*DEFF_NODE*/

static int deff_node(int **n, int *nj, int llx, int lly) {

int i = 0, j = 0, nonzero = 0;

  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nonzero += (n[i][j] > 0);
  for (j = 0; j < lly; j++)
    nonzero -= (nj[j] > 0);

  return nonzero;

}/*DEFF_NODE*/

static double entropy1d(int *n, int llx, int num) {

int i = 0;
double res = 0, alpha = 1 / (double)llx;

  for (i = 0; i < llx; i++)
      res += ((double)(n[i] + 1) / (double)(num + llx) - alpha) *
               log((double)(n[i] + 1) / (double)(num + llx));

  return res;

}/*ENTROPY1D*/

static double entropy2d(int **n, int *nj, int llx, int lly, int num) {

int i = 0, j = 0;
double res = 0, alpha = 1 / (double)(llx * lly);

  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      res += ((double)(n[i][j] + 1) / ((double)(num + llx * lly)) - alpha) *
               log((double)(n[i][j] + 1) / (double)(nj[j] + llx));

  return res;

}/*ENTROPY2D*/

/* Harald Steck optimal imaginary sample size estimator. */
SEXP alpha_star(SEXP x, SEXP data, SEXP debug) {

int i = 0, nnodes = length(data), nobs = length(VECTOR_ELT(data, 0));
int **columns = NULL, *levels = NULL;
int *node_data = NULL;
bool debugging = isTRUE(debug);
long double num = 0, den = 0, ratio = 0;
counts1d marginal = { 0 };
counts2d joint = { 0 };
SEXP nodes, labels, cur, node_info, parents, par_data, temp, cfg;

  /* dereference the columns of the data frame. */
  columns = Calloc1D(nnodes, sizeof(int *));
  levels = Calloc1D(nnodes, sizeof(int));
  for (i = 0; i < nnodes; i++) {

    temp = VECTOR_ELT(data, i);
    columns[i] = INTEGER(temp);
    levels[i] = NLEVELS(temp);

  }/*FOR*/

  /* iterate over the nodes in the graph. */
  nodes = getListElement(x, "nodes");
  PROTECT(labels = getAttrib(nodes, R_NamesSymbol));

  for (i = 0; i < nnodes; i++) {

    /* get the node label. */
    PROTECT(cur = mkString(CHAR(STRING_ELT(labels, i))));

    if (debugging)
      Rprintf("* processing node %s.\n", CHAR(STRING_ELT(cur, 0)));

    /* extract the corresponding variable from the data. */
    node_data = INTEGER(c_dataframe_column(data, cur, TRUE, FALSE));
    /* create and initialize the contingency table. */
    node_info = getListElement(nodes, (char *)CHAR(STRING_ELT(cur, 0)));
    parents = getListElement(node_info, "parents");

    if (length(parents) == 0) {

      marginal = new_1d_table(levels[i]);
      fill_1d_table(node_data, &marginal, nobs);

      /* add the effective degrees of freedom to the numerator. */
      num += deff_node_root(marginal.n, marginal.llx);
      /* add the entropy delta to the denominator. */
      den += entropy1d(marginal.n, marginal.llx, marginal.nobs);

      Free1DTAB(marginal);

    }/*THEN*/
    else {

      /* compute the parents configurations. */
      PROTECT(par_data = c_dataframe_column(data, parents, FALSE, FALSE));
      PROTECT(cfg = c_configurations(par_data, TRUE, TRUE));

      joint = new_2d_table(levels[i], NLEVELS(cfg), TRUE);
      fill_2d_table(node_data, INTEGER(cfg), &joint, nobs);

      /* add the effective degrees of freedom to the numerator. */
      num += deff_node(joint.n, joint.nj, joint.llx, joint.lly);
      /* add the entropy delta to the denominator. */
      den += entropy2d(joint.n, joint.nj, joint.llx, joint.lly, joint.nobs);

      Free2DTAB(joint);

      UNPROTECT(2);

    }/*ELSE*/

    if (debugging)
      Rprintf("  > numerator is now %Lf, denominator is now %Lf.\n", num, den);

    UNPROTECT(1);

  }/*FOR*/

  UNPROTECT(1);
  Free1D(columns);
  Free1D(levels);

  /* check that the proposed alpha is equal or greater than 1. */
  ratio = num / den;

  if (ratio < 1)
    warning("the estimated imaginary sample size is very small (%Lf).", ratio);

  return ScalarReal((double)ratio);

}/*ALPHA_STAR*/
