#include "../include/rcore.h"
#include "../core/allocations.h"
#include "../core/data.table.h"
#include "../core/sampling.h"
#include "../core/sets.h"
#include "../fitted/fitted.h"
#include "../include/sampling.h"
#include "../math/hyperpoisson.h"
#include "../math/xnegbin.h"
#include "../scores/scores.h"

static void rbn_discrete_fixed(fixed_node fixed, int *gen, int num);
static void rbn_discrete_root(ldist ld, int *gen, int num, fixed_node fixed);
static void rbn_discrete_cond(ldist ld, tabular result, int *gen,
    fixed_node fixed, bool debugging);
static void rbn_numeric_fixed(fixed_node fixed, double *gen, int num);
static void rbn_gaussian(ldist ld, tabular result, double *gen,
    fixed_node fixed);
static void rbn_mixedcg(ldist ld, tabular result, double *gen,
    fixed_node fixed);
static void rbn_zihp(ldist ld, tabular result, double *gen,
    fixed_node fixed);
static void rbn_zinb(ldist ld, tabular result, double *gen,
    fixed_node fixed);

void c_rbn_master(fitted_bn fitted, tabular result, fixed_node *fixed,
    bool debugging) {

int i = 0, k = 0, cur = 0, *poset = NULL, nnodes = fitted.nnodes;
int num = result.m.nobs;
fitted_node_e cur_node_type = ENOFIT;
ldist ld;

  /* order the nodes according to their depth in the graph. */
  poset = Calloc1D(nnodes, sizeof(int));
  topological_sort_bn(fitted, poset);

  if (debugging) {

    Rprintf("* partial node ordering is:");

    for (i = 0; i < nnodes; i++)
      Rprintf(" %s", fitted.labels[poset[i]]);

    Rprintf(".\n");

  }/*THEN*/

  /* initialize the random number generator. */
  GetRNGstate();

  for (i = 0; i < nnodes; i++) {

    /* get the index of the node we have to generate random observations from,
     * its local distribution and its type. */
    cur = poset[i];
    ld = fitted.ldists[cur];
    cur_node_type = fitted.node_types[cur];

    if (debugging) {

      if (fixed[cur].fixed)
        Rprintf("* node %s is fixed%s.\n", fitted.labels[cur],
          (ld.nparents > 0) ? ", ignoring parents" : "");
      else if (ld.nparents == 0)
        Rprintf("* simulating node %s, which doesn't have any parent.\n",
          fitted.labels[cur]);
      else {

        Rprintf("* simulating node %s with parents ", fitted.labels[cur]);
        for (k = 0; k < ld.nparents - 1; k++)
          Rprintf("%s, ", fitted.labels[ld.parents[k]]);
        Rprintf("%s.\n", fitted.labels[ld.parents[ld.nparents - 1]]);

      }/*ELSE*/

    }/*THEN*/

    /* generate the random observations for the current node. */
    switch(cur_node_type) {

      case DNODE:
      case ONODE: {

        int *gen = result.dcol[result.map[cur]];

        if (ld.nparents == 0)
          rbn_discrete_root(ld, gen, num, fixed[cur]);
        else
          rbn_discrete_cond(ld, result, gen, fixed[cur], debugging);

        break;

      }/*DNODE*/

      case GNODE: {

        double *gen = result.ccol[result.map[cur]];

        rbn_gaussian(ld, result, gen, fixed[cur]);

        break;

      }/*GNODE*/

      case CGNODE: {

        double *gen = result.ccol[result.map[cur]];

        rbn_mixedcg(ld, result, gen, fixed[cur]);

        break;

      }/*CGNODE*/

      case ZIHPNODE: {

        double *gen = result.ccol[result.map[cur]];

        rbn_zihp(ld, result, gen, fixed[cur]);

        break;

      }/*ZIHPNODE*/

      case ZINBNODE: {

        double *gen = result.ccol[result.map[cur]];

        rbn_zinb(ld, result, gen, fixed[cur]);

        break;

      }/*ZINBNODE*/

      default:
        error("unknown node type.");

    }/*SWITCH*/

  }/*FOR*/

  PutRNGstate();

  Free1D(poset);

}/*C_RBN_MASTER*/

/* sampling from a fixed discrete node (evidence). */
static void rbn_discrete_fixed(fixed_node fixed, int *gen, int num) {

  if (fixed.n == 1) {

    /* conditioning on a single level. */
    for (int i = 0; i < num; i++)
      gen[i] = fixed.idx[0];

  }/*THEN*/
  else {

    /* conditioning on a set of levels, sampled uniformly. */
    SampleReplace(num, fixed.n, gen, fixed.idx);

  }/*ELSE*/

}/*RBN_DISCRETE_FIXED*/

/* unconditional discrete sampling. */
static void rbn_discrete_root(ldist ld, int *gen, int num, fixed_node fixed) {

int np = ld.d.dims[0], *workplace = NULL;
double *p = NULL;

  if (fixed.fixed) {

    rbn_discrete_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    workplace = Calloc1D(np, sizeof(int));

    /* duplicate the probability table to save the original copy from tampering. */
    p = Calloc1D(np, sizeof(double));
    memcpy(p, ld.d.cpt, np * sizeof(double));

    /* perform the random sampling. */
    ProbSampleReplace(np, p, workplace, num, gen);

    Free1D(workplace);
    Free1D(p);

  }/*ELSE*/

}/*RBN_DISCRETE_ROOT*/

/* conditional discrete sampling. */
static void rbn_discrete_cond(ldist ld, tabular result, int *gen,
    fixed_node fixed, bool debugging) {

int nlevels = ld.d.dims[0], nconfigs = ld.d.nconfigs, np = nlevels * nconfigs;
int num = result.m.nobs;
int *workplace = NULL, *config = NULL, **dcol = NULL, *nlvls = NULL;
double *p = NULL;
bool warn = FALSE;

  if (fixed.fixed) {

    rbn_discrete_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    /* gather the (discrete) parent columns and their numbers of levels. */
    dcol = Calloc1D(ld.nparents, sizeof(int *));
    nlvls = Calloc1D(ld.nparents, sizeof(int));
    for (int k = 0; k < ld.nparents; k++) {

      dcol[k] = result.dcol[result.map[ld.parents[k]]];
      nlvls[k] = result.nlvl[result.map[ld.parents[k]]];

    }/*FOR*/

    /* allocate and initialize the parents' configurations. */
    config = Calloc1D(num, sizeof(int));
    c_fast_config(dcol, num, ld.nparents, nlvls, config, NULL, 0);

    workplace = Calloc1D(np, sizeof(int));

    /* duplicate the probability table to save the original copy from tampering. */
    p = Calloc1D(np, sizeof(double));
    memcpy(p, ld.d.cpt, np * sizeof(double));
    /* perform the random sampling. */
    CondProbSampleReplace(nlevels, nconfigs, p, config, workplace, num, gen,
      &warn);

    Free1D(dcol);
    Free1D(nlvls);
    Free1D(config);
    Free1D(workplace);
    Free1D(p);

  }/*ELSE*/

  /* warn when returning missing values. */
  if (warn && debugging)
    Rprintf("  > some parents configurations have undefined conditional distributions, NAs will be generated.");

}/*RBN_DISCRETE_COND*/

/* sampling from a fixed numeric node (evidence). */
static void rbn_numeric_fixed(fixed_node fixed, double *gen, int num) {

int i = 0;

  if (fixed.n == 1) {

    /* conditioning on a single value. */
    for (i = 0; i < num; i++)
      gen[i] = fixed.val[0];

  }/*THEN*/
  else {

    double offset = fixed.val[0], range = fixed.val[1] - fixed.val[0];

    /* conditioning on an interval, picking a value at random
     * from a uniform distribution. */
    for (i = 0; i < num; i++)
      gen[i] = offset + unif_rand() * range;

  }/*ELSE*/

}/*RBN_NUMERIC_FIXED*/

/* conditional and unconditional normal sampling. */
static void rbn_gaussian(ldist ld, tabular result, double *gen,
    fixed_node fixed) {

int i = 0, j = 0, p = ld.g.ncoefs, num = result.m.nobs;
double *beta = ld.g.coefs, sd = ld.g.sd, *Xj = NULL;

  if (fixed.fixed) {

    rbn_numeric_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    /* initialize with intercept and standard error. */
    for (i = 0; i < num; i++)
      gen[i] = beta[0] + norm_rand() * sd;

    /* add the contributions of the other regressors (if any). */
    for (j = 1; j < p; j++) {

      Xj = result.ccol[result.map[ld.parents[j - 1]]];

      for (i = 0; i < num; i++)
        gen[i] += Xj[i] * beta[j];

    }/*FOR*/

  }/*ELSE*/

}/*RBN_GAUSSIAN*/

/* conditional linear Gaussian sampling. */
static void rbn_mixedcg(ldist ld, tabular result, double *gen,
    fixed_node fixed) {

int i = 0, j = 0, ndp = ld.cg.ndparents, ngp = ld.cg.ngparents;
int num = result.m.nobs;
double *beta = ld.cg.coefs, *sd = ld.cg.sd, *beta_offset = NULL;

  if (fixed.fixed) {

    rbn_numeric_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    int **dcol = NULL, *nlvls = NULL, *config = NULL, config_nlvl = 0;
    double **ccol = NULL;

    /* gather the discrete and continuous parents. */
    ccol = Calloc1D(ngp, sizeof(double *));
    dcol = Calloc1D(ndp, sizeof(int *));
    nlvls = Calloc1D(ndp, sizeof(int));

    for (i = 0; i < ngp; i++)
      ccol[i] = result.ccol[result.map[ld.cg.gparents[i]]];

    for (i = 0; i < ndp; i++) {

      dcol[i] = result.dcol[result.map[ld.cg.dparents[i]]];
      nlvls[i] = result.nlvl[result.map[ld.cg.dparents[i]]];

    }/*FOR*/

    /* generate configurations from the discrete parents. */
    config = Calloc1D(num, sizeof(int));
    c_fast_config(dcol, num, ndp, nlvls, config, &config_nlvl, 0);

    for (i = 0; i < num; i++) {

      /* if the configuration is missing, the random observation is also a
       * missing value. */
      if (config[i] == NA_INTEGER) {

        gen[i] = NA_REAL;
        continue;

      }/*THEN*/

      /* get the right set of coefficients based on the configuration of the
       * discrete parents. */
      beta_offset = beta + (ngp + 1) * config[i];
      /* initialize with intercept and standard error. */
      gen[i] = beta_offset[0] + norm_rand() * sd[config[i]];

      for (j = 0; j < ngp; j++)
        gen[i] += ccol[j][i] * beta_offset[j + 1];

    }/*FOR*/

    Free1D(ccol);
    Free1D(dcol);
    Free1D(nlvls);
    Free1D(config);

  }/*ELSE*/

}/*RBN_MIXEDCG*/

/* build a borrowing tabular view of the continuous parents of a node. */
static tabular rbn_parents_tabular(ldist ld, tabular result) {

tabular data = empty_tabular(result.m.nobs, 0, ld.nparents);

  for (int k = 0; k < ld.nparents; k++) {

    data.ccol[k] = result.ccol[result.map[ld.parents[k]]];
    data.map[k] = k;
    data.m.flag[k].continuous = TRUE;

  }/*FOR*/

  return data;

}/*RBN_PARENTS_TABULAR*/

/* conditional and unconditional zero-inflated Poisson sampling. */
static void rbn_zihp(ldist ld, tabular result, double *gen,
    fixed_node fixed) {

int i = 0, num = result.m.nobs;
double *zinf_coefs = ld.zihp.inflation, *inten_coefs = ld.zihp.intensity;
double disp = ld.zihp.dispersion;

  if (fixed.fixed) {

    rbn_numeric_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    double *zinf_prob = Calloc1D(num, sizeof(double));
    double *lambda = Calloc1D(num, sizeof(double));
    tabular data = rbn_parents_tabular(ld, result);

    /* translate the regression coefficients into canonical parameters. */
    zihp_coefs_to_params(zinf_prob, zinf_coefs, lambda, inten_coefs, data);

    for (i = 0; i < num; i++) {

      /* set a structural zero with the prescribed probability. */
      if (unif_rand() < zinf_prob[i])
        gen[i] = 0;
      else
        gen[i] = rhypois(lambda[i], disp);

    }/*FOR*/

    FreeTAB(data);
    Free1D(zinf_prob);
    Free1D(lambda);

  }/*ELSE*/

}/*RBN_ZIHP*/

/* conditional and unconditional zero-inflated negative binomial sampling. */
static void rbn_zinb(ldist ld, tabular result, double *gen,
    fixed_node fixed) {

int i = 0, num = result.m.nobs;
double *zinf_coefs = ld.zinb.inflation, *succ_coefs = ld.zinb.prsucc;
double fail_num = ld.zinb.failures;

  if (fixed.fixed) {

    rbn_numeric_fixed(fixed, gen, num);

  }/*THEN*/
  else {

    double *zinf_prob = Calloc1D(num, sizeof(double));
    double *succ_prob = Calloc1D(num, sizeof(double));
    tabular data = rbn_parents_tabular(ld, result);

    /* translate the regression coefficients into canonical parameters. */
    zinb_coefs_to_params(zinf_prob, zinf_coefs, succ_prob, succ_coefs, data);

    for (i = 0; i < num; i++) {

      /* set a structural zero with the prescribed probability. */
      if (unif_rand() < zinf_prob[i])
        gen[i] = 0;
      else
        gen[i] = rxnegbin(succ_prob[i], fail_num);

    }/*FOR*/

    FreeTAB(data);
    Free1D(zinf_prob);
    Free1D(succ_prob);

  }/*ELSE*/

}/*RBN_ZINB*/
