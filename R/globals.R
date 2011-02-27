
# Global variables.
available.discrete.tests = c("mi", "aict", "mi-sh", "x2", "mc-mi", "mc-x2")
available.continuous.tests = c("cor", "zf", "mi-g", "mi-g-sh", "mc-mi-g", "mc-cor", "mc-zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

resampling.tests = c("mc-mi", "mc-x2", "mc-mi-g", "mc-cor", "mc-zf")
asymptotic.tests = c("mi", "mi-g", "x2", "zf")

available.discrete.scores = c("loglik", "aic", "bic", "bde", "k2")
available.continuous.scores = c("bge", "loglik-g", "aic-g", "bic-g")
available.scores = c(available.discrete.scores, available.continuous.scores)

score.equivalent.scores = c("loglik", "aic", "bic", "bde", "bge", "loglik-g", "aic-g", "bic-g")

markov.blanket.algorithms = c("gs", "iamb", "fast.iamb", "inter.iamb")
local.search.algorithms = c("mmpc")
constraint.based.algorithms = c(markov.blanket.algorithms, local.search.algorithms)
score.based.algorithms = c("hc", "tabu")
hybrid.algorithms = c("rsmax2", "mmhc")
available.learning.algorithms = c(constraint.based.algorithms, score.based.algorithms, hybrid.algorithms)

always.dag.result = c(score.based.algorithms, hybrid.algorithms)

available.mvber.vartests = c("tvar", "gvar", "nvar", "nvark")

method.labels = c(
  'gs' = "Grow-Shrink",
  'iamb' = "Incremental Association",
  'fast.iamb' = "Fast Incremental Association",
  'inter.iamb' = "Interleaved Incremental Association",
  'rnd' = "random/generated",
  'hc' = "Hill-Climbing",
  'tabu' = "Tabu Search",
  'mmpc' = "Max-Min Parent Children",
  'rsmax2' = "Two-Phase Restricted Maximization",
  'mmhc' = "Max-Min Hill-Climbing"
)

method.extra.args = list(
  'hc' = c("max.iter", "restart", "perturb"),
  'tabu' = c("max.iter", "tabu", "max.tabu")
)

test.labels = c(
  'mi' = "Mutual Information (discrete)",
  'mi-sh' = "Mutual Information (discrete, shrinkage)",
  'mc-mi' = "Mutual Information (discrete, Monte Carlo)",
  'mi-g' = "Mutual Information (Gaussian)",
  'mi-g-sh' = "Mutual Information (Gaussian, shrinkage)",
  'mc-mi-g' = "Mutual Information (Gaussian, Monte Carlo)",
  'aict'= "AIC-like Test",
  'x2'= "Pearson's X^2",
  'mc-x2'= "Pearson's X^2 (Monte Carlo)",
  'cor' = "Pearson's Linear Correlation",
  'mc-cor' = "Pearson's Linear Correlation (Monte Carlo)",
  'zf' = "Fisher's Z Test",
  'mc-zf' = "Fisher's Z Test (Monte Carlo)"
)

score.labels = c(
  'k2' = "Cooper & Herskovits' K2",
  'bde' = "Bayesian Dirichlet (BDeu)",
  'aic' = "Akaike Information Criterion",
  'bic' = "Bayesian Information Criterion",
  'loglik' = "Log-Likelihood",
  'bge' = "Bayesian Gaussian (BGe)",
  'loglik-g' = "Log-Likelihood",
  'aic-g' = "Akaike Information Criterion (Gaussian)",
  'bic-g' = "Bayesian Information Criterion (Gaussian)"
)

score.extra.args = list(
  "k2" = character(0),
  "bde" = "iss",
  "aic" = "k",
  "bic" = "k",
  "bge" = c("iss", "phi"),
  "loglik" = character(0),
  "loglik-g" = character(0),
  "aic" = "k",
  "bic" = "k",
  "aic-g" = "k",
  "bic-g" = "k"

)

graph.generation.algorithms = c("ordered", "ic-dag", "melancon", "empty", "naive")

graph.generation.labels = c(
  "ordered" = "Full Ordering",
  "ic-dag" = "Ide & Cozman's Multiconnected DAGs",
  "melancon" = "Melancon's Uniform Probability DAGs",
  "empty" = "Empty",
  "naive" = "Naive Bayes Classifier"
)

graph.generation.extra.args = list(
  "ordered" = "prob",
  "ic-dag" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree"),
  "melancon" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree")
)

cpq.algorithms = c("ls")

cpq.labels = c(
  "ls" = "Logic/Forward Sampling"
)

cpq.extra.args = list(
  "ls" = c("n", "batch")
)

discrete.loss.functions = c("logl", "pred")
continuous.loss.functions = c("logl-g")
loss.functions = c(discrete.loss.functions, continuous.loss.functions)

loss.labels = c(
  "logl" = "Log-Likelihood Loss (discrete)",
  "pred" = "Classification Error",
  "logl-g" = "Log-Likelihood Loss (Gaussian)"
)

loss.extra.args = list(
  "logl" = character(0),
  "pred" = "target",
  "logl-g" = character(0)
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

mvber.labels = list(
  "tvar" = "Total Variance",
  "gvar" = "Generalized Variance",
  "nvar" = "Squared Frobenius Norm (1/4)",
  "nvark" = "Squared Frobenius Norm (k/4)"
)

graphviz.enabled = FALSE
lattice.enabled = FALSE

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
  "hartemink" = c("initial.breaks")
)
