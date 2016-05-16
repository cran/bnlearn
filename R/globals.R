
# Global variables.
available.discrete.tests = c("mi", "mi-sh", "x2", "mc-mi", "smc-mi", "mi-adf",
  "x2-adf", "mc-x2", "smc-x2", "sp-mi", "sp-x2")
available.ordinal.tests = c("jt", "mc-jt", "smc-jt")
available.continuous.tests = c("cor", "zf", "mi-g", "mi-g-sh", "mc-mi-g",
  "smc-mi-g", "mc-cor", "smc-cor", "mc-zf", "smc-zf")
available.mixedcg.tests = c("mi-cg")
available.tests = c(available.discrete.tests, available.ordinal.tests,
  available.continuous.tests, available.mixedcg.tests)

semiparametric.tests = c("sp-mi", "sp-x2")
resampling.tests = c("mc-mi", "smc-mi", "mc-x2", "smc-x2", "mc-mi-g", "smc-mi-g",
  "mc-cor", "smc-cor", "mc-zf", "smc-zf", "mc-jt", "smc-jt", semiparametric.tests)
asymptotic.tests = c("mi", "mi-adf", "mi-g", "x2", "x2-adf", "zf", "jt", "mi-sh",
  "mi-g-sh")

available.discrete.scores = c("loglik", "aic", "bic", "bde", "bds", "k2", "mbde")
available.continuous.scores = c("loglik-g", "aic-g", "bic-g", "bge")
available.mixedcg.scores = c("loglik-cg", "aic-cg", "bic-cg")
available.scores = c(available.discrete.scores, available.continuous.scores,
  available.mixedcg.scores)

available.discrete.mi = c("mi")
available.continuous.mi = c("mi-g")
available.mi = c(available.discrete.mi, available.continuous.mi)

markov.blanket.algorithms = c("gs", "iamb", "fast.iamb", "inter.iamb")
local.search.algorithms = c("mmpc", "si.hiton.pc")
constraint.based.algorithms = c(markov.blanket.algorithms, local.search.algorithms)
score.based.algorithms = c("hc", "tabu")
hybrid.algorithms = c("rsmax2", "mmhc")
mim.based.algorithms = c("chow.liu", "aracne")
classifiers = c("naive.bayes", "tree.bayes")
available.learning.algorithms = c(constraint.based.algorithms, score.based.algorithms,
  hybrid.algorithms, mim.based.algorithms, classifiers)

always.dag.result = c(score.based.algorithms, hybrid.algorithms, classifiers)

method.labels = c(
  'gs' = "Grow-Shrink",
  'iamb' = "IAMB",
  'fast.iamb' = "Fast-IAMB",
  'inter.iamb' = "Inter-IAMB",
  'rnd' = "random/generated",
  'hc' = "Hill-Climbing",
  'tabu' = "Tabu Search",
  'mmpc' = "Max-Min Parent Children",
  'si.hiton.pc' = "Semi-Interleaved HITON-PC",
  'rsmax2' = "Two-Phase Restricted Maximization",
  'mmhc' = "Max-Min Hill-Climbing",
  'aracne' = "ARACNE",
  'chow.liu' = "Chow-Liu",
  "naive.bayes" = "Naive Bayes Classifier",
  "tree.bayes"   = "TAN Bayes Classifier"
)

method.extra.args = list(
  'hc' = c("max.iter", "maxp", "restart", "perturb"),
  'tabu' = c("max.iter", "maxp", "tabu", "max.tabu")
)

test.labels = c(
  'mi' = "Mutual Information (disc.)",
  'mi-adf' = "Mutual Information (disc., adj. d.f.)",
  'mi-sh' = "Mutual Information (disc., shrink.)",
  'mc-mi' = "Mutual Information (disc., MC)",
  'smc-mi' = "Mutual Information (disc., Seq. MC)",
  'sp-mi' = "Mutual Information (disc., semipar.)",
  'mi-g' = "Mutual Information (Gauss.)",
  'mi-g-sh' = "Mutual Information (Gauss., shrink.)",
  'mc-mi-g' = "Mutual Information (Gauss., MC)",
  'smc-mi-g' = "Mutual Information (Gauss., Seq. MC)",
  'mi-cg' = "Mutual Information (cond. Gauss.)",
  'x2' = "Pearson's X^2",
  'x2-adf' = "Pearson's X^2 (adj. d.f.)",
  'mc-x2'= "Pearson's X^2 (MC)",
  'smc-x2'= "Pearson's X^2 (Seq. MC)",
  'sp-x2'= "Pearson's X^2 (semipar.)",
  'jt' = "Jonckheere-Terpstra",
  'mc-jt' = "Jonckheere-Terpstra (MC)",
  'smc-jt' = "Jonckheere-Terpstra (Seq. MC)",
  'cor' = "Pearson's Correlation",
  'mc-cor' = "Pearson's Correlation (MC)",
  'smc-cor' = "Pearson's Correlation (Seq. MC)",
  'zf' = "Fisher's Z",
  'mc-zf' = "Fisher's Z (MC)",
  'smc-zf' = "Fisher's Z (Seq. MC)"
)

score.labels = c(
  'k2' = "Cooper & Herskovits' K2",
  'bde' = "Bayesian Dirichlet (BDe)",
  'bds' = "Bayesian Dirichlet Sparse (BDs)",
  'mbde' = "Bayesian Dirichlet (interventional data)",
  'aic' = "AIC (disc.)",
  'bic' = "BIC (disc.)",
  'loglik' = "Log-Likelihood (disc.)",
  'bge' = "Bayesian Gaussian (BGe)",
  'loglik-g' = "Log-Likelihood (Gauss.)",
  'aic-g' = "AIC (Gauss.)",
  'bic-g' = "BIC (Gauss.)",
  'loglik-cg' = "Log-Likelihood (cond. Gauss.)",
  'aic-cg' = "AIC (cond. Gauss.)",
  'bic-cg' = "BIC (cond. Gauss.)"
)

score.extra.args = list(
  "k2" = character(0),
  "bde" = c("prior", "beta", "iss"),
  "bds" = c("prior", "beta", "iss"),
  "mbde" = c("prior", "beta", "iss", "exp"),
  "aic" = c("k"),
  "bic" = c("k"),
  "bge" = c("prior", "beta", "iss", "phi"),
  "loglik" = character(0),
  "loglik-g" = character(0),
  "aic-g" = c("k"),
  "bic-g" = c("k"),
  "loglik-cg" = character(0),
  "aic-cg" = c("k"),
  "bic-cg" = c("k")
)

mi.estimator.labels = c(
  'mi' = "Maximum Likelihood (disc.)",
  'mi-g' = "Maximum Likelihood (Gauss.)"
)

mi.estimator.tests = c(
  'mi' = "mi",
  'mi-g' = "mi-g"
)

graph.generation.algorithms = c("ordered", "ic-dag", "melancon", "empty", "averaged")

graph.generation.labels = c(
  "ordered" = "Full Ordering",
  "ic-dag" = "Ide & Cozman's Multiconnected DAGs",
  "melancon" = "Melancon's Uniform Probability DAGs",
  "empty" = "Empty",
  "averaged" = "Model Averaging"
)

graph.generation.extra.args = list(
  "ordered" = "prob",
  "ic-dag" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree", "every"),
  "melancon" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree", "every"),
  "averaged" = "threshold"
)

prior.distributions = c("uniform", "vsp", "cs")

prior.labels = c(
  "uniform" = "Uniform",
  "vsp" = "Variable Selection",
  "cs" = "Castelo & Siebes"
)

cpq.algorithms = c("ls", "lw")

cpq.labels = c(
  "ls" = "Logic/Forward Sampling",
  "lw" = "Likelihood Weighting"
)

cpq.extra.args = list(
  "ls" = c("n", "batch", "query.nodes"),
  "lw" = c("n", "batch", "query.nodes")
)

discrete.loss.functions = c("logl", "pred", "pred-lw")
continuous.loss.functions = c("logl-g", "cor", "cor-lw", "mse", "mse-lw")
mixedcg.loss.functions = c("logl-cg", "cor-lw-cg", "mse-lw-cg", "pred-lw-cg")
loss.functions = c(discrete.loss.functions, continuous.loss.functions,
  mixedcg.loss.functions)

loss.labels = c(
  "logl" = "Log-Likelihood Loss (disc.)",
  "pred" = "Classification Error",
  "pred-lw" = "Classification Error (Posterior, disc.)",
  "pred-lw-cg" = "Classification Error (Posterior, cond. Gauss.)",
  "logl-g" = "Log-Likelihood Loss (Gauss.)",
  "cor" = "Predictive Correlation",
  "cor-lw" = "Predictive Correlation (Posterior, Gauss.)",
  "cor-lw-cg" = "Predictive Correlation (Posterior, cond. Gauss.)",
  "mse" = "Mean Squared Error",
  "mse-lw" = "Mean Squared Error (Posterior, Gauss.)",
  "mse-lw-cg" = "Mean Squared Error (Posterior, cond. Gauss.)",
  "logl-cg" = "Log-Likelihood Loss (cond. Gauss.)"
)

loss.extra.args = list(
  "logl" = character(0),
  "pred" = "target",
  "pred-lw" = c("target", "n", "from"),
  "pred-lw-cg" = c("target", "n", "from"),
  "logl-g" = character(0),
  "cor" = "target",
  "cor-lw" = c("target", "n", "from"),
  "cor-lw-cg" = c("target", "n", "from"),
  "mse" = "target",
  "mse-lw" = c("target", "n", "from"),
  "mse-lw-cg" = c("target", "n", "from"),
  "logl-cg" = character(0)
)

available.fitting.methods = c("mle", "bayes")

fitting.labels = c(
  "mle" = "Maximum Likelihood",
  "bayes" = "Bayesian Parameter Estimation"
)

fitting.extra.args = list(
  "mle" = character(0),
  "bayes" = "iss"
)

available.prediction.methods = c("parents", "bayes-lw")

prediction.labels = c(
  "parents" = "Parents (Maximum Likelihood)",
  "bayes-lw" = "Posterior Expectation (Likelihood Weighting)"
)

prediction.extra.args = list(
  "parents" = character(0),
  "bayes-lw" = c("n", "from")
)

supported.clusters = c("MPIcluster", "PVMcluster","SOCKcluster")

available.discretization.methods = c("quantile", "interval", "hartemink")

discretization.labels = c(
  "quantile" = "Quantile Discretization",
  "interval" = "Interval Discretization",
  "hartemink" = "Hartemink's Pairwise Mutual Information"
)

discretization.extra.args = list(
  "quantile" = character(0),
  "interval" = character(0),
  "hartemink" = c("ibreaks", "idisc")
)

fitted.from.data = c(
  "continuous" = "bn.fit.gnet",
  "factor" = "bn.fit.dnet",
  "ordered" = "bn.fit.onet",
  "mixed-cg" = "bn.fit.cgnet",
  "mixed-do" = "bn.fit.donet"
)

discrete.data.types = c("factor", "ordered", "mixed-do")
continuous.data.types = c("continuous")
mixed.data.types = c("mixed-cg")
available.data.types = c(discrete.data.types, continuous.data.types,
  mixed.data.types)

data.type.labels = c(
  "continuous" = "all variables must be numeric",
  "factor" = "all variables must be unordered factors",
  "ordered" = "all variables must be ordered factors",
  "mixed-cg" = "variables can be either numeric or factors",
  "mixed-do" = "variables can be either ordered or unordered factors"
)

fitted.node.types = c("bn.fit.dnode", "bn.fit.onode", "bn.fit.gnode",
  "bn.fit.cgnode")

# global test counter.
reset.test.counter = function() {

  invisible(.Call("reset_test_counter"))

}#RESET.TEST.COUNTER

increment.test.counter = function(i = 1) {

  if (!is.real.number(i))
    stop("the increment must be a single real number.")

  invisible(.Call("increment_test_counter", i))

}#INCREMENT.TEST.COUNTER

test.counter = function() {

  return(.Call("get_test_counter"))

}#TEST.COUNTER

