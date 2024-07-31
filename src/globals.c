#include "include/rcore.h"
#include "include/register.h"
#include <R_ext/Rdynload.h>

SEXP BN_ModelstringSymbol;
SEXP BN_NodesSymbol;
SEXP BN_ProbSymbol;
SEXP BN_MethodSymbol;
SEXP BN_WeightsSymbol;
SEXP BN_DsepsetSymbol;
SEXP BN_MetaDataSymbol;
SEXP BN_NobsSymbol;
SEXP BN_DfSymbol;
SEXP TRUESEXP, FALSESEXP;

double *regret_table;

SEXP onLoad(void) {

  /* initialize symbols in .onLoad(), to do that only once. */
  BN_ModelstringSymbol = install("modelstring");
  BN_NodesSymbol = install("nodes");
  BN_ProbSymbol = install("prob");
  BN_MethodSymbol = install("method");
  BN_WeightsSymbol = install("weights");
  BN_DsepsetSymbol = install("dsep.set");
  BN_MetaDataSymbol = install("metadata");
  BN_NobsSymbol = install("nobs");
  BN_DfSymbol = install("df");
  TRUESEXP = ScalarLogical(TRUE);
  R_PreserveObject(TRUESEXP);
  FALSESEXP = ScalarLogical(FALSE);
  R_PreserveObject(FALSESEXP);

  return R_NilValue;

}/*C_ONLOAD*/

SEXP onUnload(void) {

  R_ReleaseObject(TRUESEXP);
  R_ReleaseObject(FALSESEXP);

  return R_NilValue;

}/*C_ONUNLOAD*/

#define CALL_ENTRY(fun, args) \
  {"call_"#fun,                (DL_FUNC) &fun,                 args}

static const R_CallMethodDef CallEntries[] = {
  CALL_ENTRY(all_equal_bn, 2),
  CALL_ENTRY(allsubs_test, 12),
  CALL_ENTRY(alpha_star, 3),
  CALL_ENTRY(amat2arcs, 2),
  CALL_ENTRY(aracne, 6),
  CALL_ENTRY(arcs_cg_assumptions, 3),
  CALL_ENTRY(arcs_rbind, 3),
  CALL_ENTRY(arcs2amat, 2),
  CALL_ENTRY(arcs2elist, 6),
  CALL_ENTRY(bn_recovery, 4),
  CALL_ENTRY(bootstrap_arc_coefficients, 2),
  CALL_ENTRY(bootstrap_reduce, 1),
  CALL_ENTRY(bootstrap_strength_counters, 4),
  CALL_ENTRY(cache_partial_structure, 4),
  CALL_ENTRY(cache_structure, 3),
  CALL_ENTRY(castelo_completion, 3),
  CALL_ENTRY(ccgpred, 4),
  CALL_ENTRY(cdpred, 4),
  CALL_ENTRY(cg_banned_arcs, 2),
  CALL_ENTRY(cgpred, 3),
  CALL_ENTRY(cgsd, 3),
  CALL_ENTRY(check_covariance, 1),
  CALL_ENTRY(chow_liu, 8),
  CALL_ENTRY(class_err, 2),
  CALL_ENTRY(classic_discrete_parameters, 6),
  CALL_ENTRY(colliders, 6),
  CALL_ENTRY(configurations, 3),
  CALL_ENTRY(count_observed_values, 1),
  CALL_ENTRY(cpdag, 9),
  CALL_ENTRY(cpdist_lw, 5),
  CALL_ENTRY(dag2ug, 3),
  CALL_ENTRY(data_frame_finite, 1),
  CALL_ENTRY(data_type, 1),
  CALL_ENTRY(dataframe_column, 4),
  CALL_ENTRY(dedup, 4),
  CALL_ENTRY(dpred, 4),
  CALL_ENTRY(elist2arcs, 1),
  CALL_ENTRY(empty_graph, 2),
  CALL_ENTRY(fit2arcs, 1),
  CALL_ENTRY(fitted_mb, 2),
  CALL_ENTRY(fitted_vs_data, 3),
  CALL_ENTRY(gaussian_ols_parameters, 6),
  CALL_ENTRY(get_test_counter, 0),
  CALL_ENTRY(gpred, 3),
  CALL_ENTRY(has_pdag_path, 8),
  CALL_ENTRY(hc_opt_step, 10),
  CALL_ENTRY(hc_to_be_added, 7),
  CALL_ENTRY(hierarchical_dirichlet_parameters, 8),
  CALL_ENTRY(ide_cozman_graph, 8),
  CALL_ENTRY(increment_test_counter, 1),
  CALL_ENTRY(indep_test, 9),
  CALL_ENTRY(is_dag, 2),
  CALL_ENTRY(is_listed, 5),
  CALL_ENTRY(is_pdag_acyclic, 5),
  CALL_ENTRY(is_row_equal, 2),
  CALL_ENTRY(joint_discretize, 7),
  CALL_ENTRY(loglikelihood_function, 7),
  CALL_ENTRY(lw_weights, 4),
  CALL_ENTRY(mappred, 7),
  CALL_ENTRY(marginal_discretize, 5),
  CALL_ENTRY(match_brace, 4),
  CALL_ENTRY(mean_strength, 3),
  CALL_ENTRY(minimal_data_frame, 1),
  CALL_ENTRY(minimal_table, 2),
  CALL_ENTRY(mixture_gaussian_ols_parameters, 7),
  CALL_ENTRY(naivepred, 7),
  CALL_ENTRY(nbr2arcs, 1),
  CALL_ENTRY(normalize_cpt, 1),
  CALL_ENTRY(nparams_cgnet, 3),
  CALL_ENTRY(nparams_fitted, 3),
  CALL_ENTRY(num_arcs, 1),
  CALL_ENTRY(onLoad, 0),
  CALL_ENTRY(onUnload, 0),
  CALL_ENTRY(ordered_graph, 3),
  CALL_ENTRY(pdag_extension, 3),
  CALL_ENTRY(pdag2dag, 2),
  CALL_ENTRY(per_node_score, 6),
  CALL_ENTRY(rbn_master, 5),
  CALL_ENTRY(reset_test_counter, 0),
  CALL_ENTRY(root_nodes, 2),
  CALL_ENTRY(roundrobin_test, 9),
  CALL_ENTRY(score_cache_fill, 13),
  CALL_ENTRY(score_delta, 9),
  CALL_ENTRY(shd, 3),
  CALL_ENTRY(smart_network_averaging, 3),
  CALL_ENTRY(subsets, 2),
  CALL_ENTRY(tabu_hash, 4),
  CALL_ENTRY(tabu_step, 13),
  CALL_ENTRY(tiers, 2),
  CALL_ENTRY(topological_ordering, 4),
  CALL_ENTRY(tree_directions, 4),
  CALL_ENTRY(unique_arcs, 3),
  CALL_ENTRY(which_undirected, 2),
  {NULL, NULL, 0}
};

void R_init_bnlearn(DllInfo *dll) {

  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

}/*R_INIT_BNLEARN*/

