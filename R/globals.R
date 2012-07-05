
# Global variables.
available.discrete.tests = c("mi", "mi-sh", "x2", "mc-mi", "smc-mi", "mc-x2", "smc-x2")
available.continuous.tests = c("cor", "zf", "mi-g", "mi-g-sh", "mc-mi-g",
  "smc-mi-g", "mc-cor", "smc-cor", "mc-zf", "smc-zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

resampling.tests = c("mc-mi", "smc-mi", "mc-x2", "smc-x2", "mc-mi-g", "smc-mi-g",
  "mc-cor", "smc-cor", "mc-zf", "smc-zf")
asymptotic.tests = c("mi", "mi-g", "x2", "zf")

available.discrete.scores = c("loglik", "aic", "bic", "bde", "bdes", "k2", "mbde")
available.continuous.scores = c("bge", "loglik-g", "aic-g", "bic-g")
available.scores = c(available.discrete.scores, available.continuous.scores)

score.equivalent.scores = c("loglik", "aic", "bic", "bde", "bge", "loglik-g",
  "aic-g", "bic-g")

available.discrete.mi = c("mi")
available.continuous.mi = c("mi-g")
available.mi = c(available.discrete.mi, available.continuous.mi)

markov.blanket.algorithms = c("gs", "iamb", "fast.iamb", "inter.iamb")
local.search.algorithms = c("mmpc", "si.hiton.pc")
constraint.based.algorithms = c(markov.blanket.algorithms, local.search.algorithms)
score.based.algorithms = c("hc", "tabu")
hybrid.algorithms = c("rsmax2", "mmhc")
mim.based.algorithms = c("chow.liu", "aracne")
classifiers = c("naive", "tan")
available.learning.algorithms = c(constraint.based.algorithms, score.based.algorithms,
  hybrid.algorithms, mim.based.algorithms, classifiers)

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
  'si.hiton.pc' = "Semi-Interleaved HITON-PC",
  'rsmax2' = "Two-Phase Restricted Maximization",
  'mmhc' = "Max-Min Hill-Climbing",
  'aracne' = "ARACNE",
  'chow.liu' = "Chow-Liu",
  "naive" = "Naive Bayes Classifier",
  "tan"   = "Tree-Augmented Naive Bayes Classifier"
)

method.extra.args = list(
  'hc' = c("max.iter", "restart", "perturb"),
  'tabu' = c("max.iter", "tabu", "max.tabu")
)

test.labels = c(
  'mi' = "Mutual Information (discrete)",
  'mi-sh' = "Mutual Information (discrete, shrinkage)",
  'mc-mi' = "Mutual Information (discrete, Monte Carlo)",
  'smc-mi' = "Mutual Information (discrete, Sequential Monte Carlo)",
  'mi-g' = "Mutual Information (Gaussian)",
  'mi-g-sh' = "Mutual Information (Gaussian, shrinkage)",
  'mc-mi-g' = "Mutual Information (Gaussian, Monte Carlo)",
  'smc-mi-g' = "Mutual Information (Gaussian, Sequential Monte Carlo)",
  'x2'= "Pearson's X^2",
  'mc-x2'= "Pearson's X^2 (Monte Carlo)",
  'smc-x2'= "Pearson's X^2 (Sequential Monte Carlo)",
  'cor' = "Pearson's Linear Correlation",
  'mc-cor' = "Pearson's Linear Correlation (Monte Carlo)",
  'smc-cor' = "Pearson's Linear Correlation (Sequential Monte Carlo)",
  'zf' = "Fisher's Z Test",
  'mc-zf' = "Fisher's Z Test (Monte Carlo)",
  'smc-zf' = "Fisher's Z Test (Sequential Monte Carlo)"
)

score.labels = c(
  'k2' = "Cooper & Herskovits' K2",
  'bde' = "Bayesian Dirichlet (BDeu)",
  'bdes' = "Sparse Bayesian Dirichlet (BDes)",
  'mbde' = "Bayesian Dirichlet (interventional data)",
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
  "bdes" = "iss",
  "mbde" = c("iss", "exp"),
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

mi.estimator.labels = c(
  'mi' = "Maximum Likelihood (discrete)",
  'mi-g' = "Maximum Likelihood (Gaussian)"
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

template.numeric = numeric(1)

