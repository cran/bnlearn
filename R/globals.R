
#-- conditional independence tests --------------------------------------------#
available.discrete.tests = c("mi", "mi-sh", "x2", "mc-mi", "smc-mi", "mi-adf",
  "x2-adf", "mc-x2", "smc-x2", "sp-mi", "sp-x2")
available.ordinal.tests = c("jt", "mc-jt", "smc-jt")
available.continuous.tests = c("cor", "zf", "mi-g", "mi-g-sh", "mc-mi-g",
  "smc-mi-g", "mc-cor", "smc-cor", "mc-zf", "smc-zf")
available.mixedcg.tests = c("mi-cg")
available.omnibus.tests = c("custom-test")
available.tests = c(available.discrete.tests, available.ordinal.tests,
  available.continuous.tests, available.mixedcg.tests, available.omnibus.tests)

semiparametric.tests = c("sp-mi", "sp-x2")
resampling.tests = c("mc-mi", "smc-mi", "mc-x2", "smc-x2", "mc-mi-g",
  "smc-mi-g", "mc-cor", "smc-cor", "mc-zf", "smc-zf", "mc-jt", "smc-jt",
  semiparametric.tests)
asymptotic.tests = c("mi", "mi-adf", "mi-g", "x2", "x2-adf", "zf", "jt",
  "mi-sh", "mi-g-sh", "mi-cg")

test.labels = c(
  "mi" = "Mutual Information (disc.)",
  "mi-adf" = "Mutual Information (disc., adj. d.f.)",
  "mi-sh" = "Mutual Information (disc., shrink.)",
  "mc-mi" = "Mutual Information (disc., MC)",
  "smc-mi" = "Mutual Information (disc., Seq. MC)",
  "sp-mi" = "Mutual Information (disc., semipar.)",
  "mi-g" = "Mutual Information (Gauss.)",
  "mi-g-sh" = "Mutual Information (Gauss., shrink.)",
  "mc-mi-g" = "Mutual Information (Gauss., MC)",
  "smc-mi-g" = "Mutual Information (Gauss., Seq. MC)",
  "mi-cg" = "Mutual Information (cond. Gauss.)",
  "x2" = "Pearson's X^2",
  "x2-adf" = "Pearson's X^2 (adj. d.f.)",
  "mc-x2"= "Pearson's X^2 (MC)",
  "smc-x2"= "Pearson's X^2 (Seq. MC)",
  "sp-x2"= "Pearson's X^2 (semipar.)",
  "jt" = "Jonckheere-Terpstra",
  "mc-jt" = "Jonckheere-Terpstra (MC)",
  "smc-jt" = "Jonckheere-Terpstra (Seq. MC)",
  "cor" = "Pearson's Correlation",
  "mc-cor" = "Pearson's Correlation (MC)",
  "smc-cor" = "Pearson's Correlation (Seq. MC)",
  "zf" = "Fisher's Z",
  "mc-zf" = "Fisher's Z (MC)",
  "smc-zf" = "Fisher's Z (Seq. MC)",
  "custom-test" = "User-Provided Test Function"
)

test.extra.args = list(
  "mi" = character(0),
  "mi-adf" = character(0),
  "mi-sh" = character(0),
  "mc-mi" = c("B"),
  "smc-mi" = c("B"),
  "sp-mi" = c("B"),
  "mi-g" = character(0),
  "mi-g-sh" = character(0),
  "mc-mi-g" = c("B"),
  "smc-mi-g" = c("B"),
  "mi-cg" = character(0),
  "x2" = character(0),
  "x2-adf" = character(0),
  "mc-x2"= c("B"),
  "smc-x2"= c("B"),
  "sp-x2"= c("B"),
  "jt" = character(0),
  "mc-jt" = c("B"),
  "smc-jt" = c("B"),
  "cor" = character(0),
  "mc-cor" = c("B"),
  "smc-cor" = c("B"),
  "zf" = character(0),
  "mc-zf" = c("B"),
  "smc-zf" = c("B"),
  "custom-test" = c("fun", "args")
)

#-- network scores ------------------------------------------------------------#
available.discrete.bayesian.scores =
  c("bde", "bds", "bdj", "k2", "mbde", "bdla")
available.discrete.scores =
  c("loglik", "aic", "bic", "ebic", "pred-loglik", "fnml", "qnml", "nal",
    "pnal", available.discrete.bayesian.scores)
available.continuous.bayesian.scores = c("bge")
available.continuous.scores =
  c("loglik-g", "aic-g", "bic-g", "ebic-g", "pred-loglik-g", "nal-g", "pnal-g",
    available.continuous.bayesian.scores)
available.mixedcg.scores =
  c("loglik-cg", "aic-cg", "bic-cg", "ebic-cg", "pred-loglik-cg", "nal-cg",
    "pnal-cg")
available.omnibus.scores = c("custom-score")
available.scores = c(available.discrete.scores, available.continuous.scores,
  available.mixedcg.scores, available.omnibus.scores)

scores.for.incomplete.data =
  c("custom-score", "nal", "pnal", "nal-g", "pnal-g", "nal-cg", "pnal-cg")

score.labels = c(
  "k2" = "Cooper & Herskovits' K2",
  "bde" = "Bayesian Dirichlet (BDe)",
  "bds" = "Bayesian Dirichlet Sparse (BDs)",
  "bdj" = "Bayesian Dirichlet, Jeffrey's prior",
  "mbde" = "Bayesian Dirichlet (interventional data)",
  "bdla" = "Bayesian Dirichlet, Locally Averaged",
  "aic" = "AIC (disc.)",
  "bic" = "BIC (disc.)",
  "ebic" = "eBIC (disc.)",
  "loglik" = "Log-Likelihood (disc.)",
  "pred-loglik" = "Predictive Log-Likelihood (disc.)",
  "fnml" = "Factorized Normalized Maximum Likelihood",
  "qnml" = "Quotient Normalized Maximum Likelihood",
  "nal" = "Node-Average Likelihood (disc.)",
  "pnal" = "Penalized Node-Average Likelihood (disc.)",
  "bge" = "Bayesian Gaussian (BGe)",
  "loglik-g" = "Log-Likelihood (Gauss.)",
  "pred-loglik-g" = "Predictive Log-Likelihood (Gauss.)",
  "aic-g" = "AIC (Gauss.)",
  "bic-g" = "BIC (Gauss.)",
  "ebic-g" = "eBIC (Gauss.)",
  "nal-g" = "Node-Average Likelihood (Gauss.)",
  "loglik-cg" = "Log-Likelihood (cond. Gauss.)",
  "pred-loglik-cg" = "Predictive Log-Likelihood (cond. Gauss.)",
  "pnal-g" = "Penalized Node-Average Likelihood (Gauss.)",
  "aic-cg" = "AIC (cond. Gauss.)",
  "bic-cg" = "BIC (cond. Gauss.)",
  "ebic-cg" = "eBIC (cond. Gauss.)",
  "nal-cg" = "Node-Average Likelihood (cond. Gauss.)",
  "pnal-cg" = "Penalized Node-Average Likelihood (cond. Gauss.)",
  "custom-score" = "User-Provided Score Function"
)

score.extra.args = list(
  "k2" = character(0),
  "bde" = c("prior", "beta", "iss"),
  "bds" = c("prior", "beta", "iss"),
  "bdj" = c("prior", "beta"),
  "mbde" = c("prior", "beta", "iss", "exp"),
  "bdla" = c("prior", "beta", "l"),
  "aic" = c("k"),
  "bic" = c("k"),
  "ebic" = c("k", "gamma"),
  "bge" = c("prior", "beta", "nu", "iss.mu", "iss.w"),
  "loglik" = character(0),
  "pred-loglik" = c("newdata"),
  "fnml" = character(0),
  "qnml" = character(0),
  "nal" = character(0),
  "pnal" = c("k"),
  "loglik-g" = character(0),
  "pred-loglik-g" = c("newdata"),
  "aic-g" = c("k"),
  "bic-g" = c("k"),
  "ebic-g" = c("k", "gamma"),
  "nal-g" = character(0),
  "pnal-g" = c("k"),
  "loglik-cg" = character(0),
  "pred-loglik-cg" = c("newdata"),
  "aic-cg" = c("k"),
  "bic-cg" = c("k"),
  "ebic-cg" = c("k", "gamma"),
  "nal-cg" = character(0),
  "pnal-cg" = c("k"),
  "custom-score" = c("fun", "args")
)

prior.distributions = c("uniform", "vsp", "cs", "marginal")

prior.labels = c(
  "uniform" = "Uniform",
  "vsp" = "Variable Selection",
  "cs" = "Castelo & Siebes",
  "marginal" = "Marginal Uniform"
)

#-- mutual information estimators ---------------------------------------------#
available.discrete.mi = c("mi")
available.continuous.mi = c("mi-g")
available.mi = c(available.discrete.mi, available.continuous.mi)

mi.estimator.labels = c(
  "mi" = "Maximum Likelihood (disc.)",
  "mi-g" = "Maximum Likelihood (Gauss.)"
)

mi.estimator.tests = c(
  "mi" = "mi",
  "mi-g" = "mi-g"
)

#-- structure learning algorithms ---------------------------------------------#
markov.blanket.algorithms =
  c("gs", "iamb", "fast.iamb", "inter.iamb", "iamb.fdr")
local.search.algorithms = c("pc.stable", "mmpc", "si.hiton.pc", "hpc")
constraint.based.algorithms =
  c(markov.blanket.algorithms, local.search.algorithms)
score.based.algorithms = c("hc", "tabu")
em.algorithms = c("structural.em")
hybrid.algorithms = c("rsmax2", "mmhc", "h2pc")
mim.based.algorithms = c("chow.liu", "aracne")
classification.algorithms = c("naive.bayes", "tree.bayes")
available.learning.algorithms = c(constraint.based.algorithms,
  score.based.algorithms, hybrid.algorithms, mim.based.algorithms,
  classification.algorithms, em.algorithms)

learning.labels = c(
  "pc.stable" = "PC (Stable)",
  "gs" = "Grow-Shrink",
  "iamb" = "IAMB",
  "fast.iamb" = "Fast-IAMB",
  "inter.iamb" = "Inter-IAMB",
  "iamb.fdr" = "IAMB-FDR",
  "rnd" = "random/generated",
  "hc" = "Hill-Climbing",
  "tabu" = "Tabu Search",
  "structural.em" = "Structural EM",
  "mmpc" = "Max-Min Parent Children",
  "si.hiton.pc" = "Semi-Interleaved HITON-PC",
  "hpc" = "Hybrid Parents and Children",
  "rsmax2" = "Two-Phase Restricted Maximization",
  "mmhc" = "Max-Min Hill-Climbing",
  "h2pc" = "Hybrid^2 Parent Children",
  "aracne" = "ARACNE",
  "chow.liu" = "Chow-Liu",
  "naive.bayes" = "Naive Bayes Classifier",
  "tree.bayes"   = "TAN Bayes Classifier"
)

learning.extra.args = list(
  "hc" = c("max.iter", "maxp", "restart", "perturb"),
  "tabu" = c("max.iter", "maxp", "tabu", "max.tabu"),
  "chow.liu" = character(0),
  "tree.bayes" = c("estimator", "root")
)

#-- random graph generation algorithms ----------------------------------------#
random.graph.generation.algorithms = c("ordered", "ic-dag", "melancon")
graph.generation.algorithms =
  c(random.graph.generation.algorithms, "empty", "complete", "averaged")

graph.generation.labels = c(
  "ordered" = "Full Ordering",
  "ic-dag" = "Ide & Cozman's Multiconnected DAGs",
  "melancon" = "Melancon's Uniform Probability DAGs",
  "empty" = "Empty",
  "complete" = "Complete DAGs",
  "averaged" = "Model Averaging"
)

graph.generation.extra.args = list(
  "ordered" = "prob",
  "ic-dag" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree",
               "every"),
  "melancon" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree",
                 "every"),
  "averaged" = "threshold"
)

#-- conditional probability query algorithms ----------------------------------#
cpq.algorithms = c("ls", "lw")

cpq.labels = c(
  "ls" = "Logic/Forward Sampling",
  "lw" = "Likelihood Weighting"
)

cpq.extra.args = list(
  "ls" = c("n", "batch", "query.nodes"),
  "lw" = c("n", "batch", "query.nodes")
)

#-- cross-validation loss functions -------------------------------------------#
discrete.loss.functions = c("logl", "pred", "pred-lw")
continuous.loss.functions = c("logl-g", "cor", "cor-lw", "mse", "mse-lw")
mixedcg.loss.functions = c("logl-cg", "cor-lw-cg", "mse-lw-cg", "pred-lw-cg")
classifiers.loss.functions = c("pred-exact")
loss.functions = c(discrete.loss.functions, continuous.loss.functions,
  mixedcg.loss.functions, classifiers.loss.functions)

loss.labels = c(
  "logl" = "Log-Likelihood Loss (disc.)",
  "pred-exact" = "Classification Error (Posterior, exact)",
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
  "pred" = c("predict", "target", "predict.args"),
  "pred-exact" = c("target", "prior"),
  "pred-lw" = c("target", "n", "from"),
  "pred-lw-cg" = c("target", "n", "from"),
  "logl-g" = character(0),
  "cor" = c("predict", "target", "predict.args"),
  "cor-lw" = c("target", "n", "from"),
  "cor-lw-cg" = c("target", "n", "from"),
  "mse" = c("predict", "target", "predict.args"),
  "mse-lw" = c("target", "n", "from"),
  "mse-lw-cg" = c("target", "n", "from"),
  "logl-cg" = character(0)
)

#-- parameter estimators ------------------------------------------------------#
available.dbn.fits = c("mle", "bayes", "hdir", "hard-em")
available.gbn.fits = c("mle-g", "hard-em-g")
available.cgbn.fits = c("mle-cg", "hard-em-cg")
available.fits = c(available.dbn.fits, available.gbn.fits, available.cgbn.fits)
complete.data.fits = c("mle", "bayes", "hdir", "mle-g", "mle-cg")

fits.labels = c(
  "mle" = "Maximum Likelihood (disc.)",
  "mle-g" = "Maximum Likelihood (Gauss.)",
  "mle-cg" = "Maximum Likelihood (cond. Gauss.)",
  "bayes" = "Bayesian Dirichlet",
  "hdir" = "Bayesian Hierarchical Dirichlet",
  "hard-em" = "Hard Expectation-Maximization (disc.)",
  "hard-em-g" = "Hard Expectation-Maximization (Gauss.)",
  "hard-em-cg" = "Hard Expectation-Maximization (cond. Gauss.)"
)

fits.extra.args = list(
  "mle" = "replace.unidentifiable",
  "mle-g" = "replace.unidentifiable",
  "mle-cg" = "replace.unidentifiable",
  "bayes" = "iss",
  "hdir" = c("iss", "alpha0", "group"),
  "hard-em" = c("impute", "impute.args", "fit", "fit.args",
                "loglik.threshold", "params.threshold", "max.iter",
                "newdata", "start"),
  "hard-em-g" = c("impute", "impute.args", "fit", "fit.args",
                  "loglik.threshold", "params.threshold", "max.iter",
                  "newdata", "start"),
  "hard-em-cg" = c("impute", "impute.args", "fit", "fit.args",
                   "loglik.threshold", "params.threshold", "max.iter",
                   "newdata", "start")
)

fitted.from.method = c(
  "mle" = "bn.fit.dnet",
  "mle-g" = "bn.fit.gnet",
  "mle-cg" = "bn.fit.cgnet",
  "bayes" = "bn.fit.dnet",
  "hdir" = "bn.fit.dnet",
  "hard-em" = "bn.fit.dnet",
  "hard-em-g" = "bn.fit.gnet",
  "hard-em-cg" = "bn.fit.cgnet"
)

#-- cross-validation fold schemes ---------------------------------------------#
available.cv.methods = c("k-fold", "hold-out", "custom-folds")

cv.labels = c(
  "k-fold" = "k-Fold",
  "hold-out" = "Hold-Out",
  "custom-folds" = "Custom Folds"
)

cv.extra.args = list(
  "k-fold" = c("k", "runs"),
  "hold-out" = c("k", "m", "runs"),
  "custom-folds" = c("folds")
)

#-- prediction methods --------------------------------------------------------#
available.prediction.methods = c("parents", "bayes-lw", "exact")

prediction.labels = c(
  "parents" = "Parents (Maximum Likelihood)",
  "bayes-lw" = "Posterior Expectation (Likelihood Weighting)",
  "exact" = "Exact Inference"
)

prediction.extra.args = list(
  "parents" = character(0),
  "bayes-lw" = c("n", "from"),
  "exact" = "from"
)

#-- imputation methods --------------------------------------------------------#
available.imputation.methods = c("parents", "bayes-lw", "exact")

imputation.extra.args = list(
  "parents" = character(0),
  "bayes-lw" = "n",
  "exact" = character(0)
)

imputation.labels = c(
  "parents" = "Parents (Maximum Likelihood)",
  "bayes-lw" = "Posterior Expectation (Likelihood Weighting)",
  "exact" = "Exact Inference"
)

#-- discretization methods ----------------------------------------------------#
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

#-- arc strength estimators ---------------------------------------------------#
available.strength.methods = c("test", "score", "bootstrap", "bayes-factor")

#-- graph enumeration formulas ------------------------------------------------#
available.enumerations = c("all-dags", "dags-disregarding-one-arc",
  "dags-given-ordering", "dags-with-k-roots", "dags-with-r-arcs")

enumerations.extra.args = list(
  "all-dags" = character(0),
  "dags-disregarding-one-arc" = character(0),
  "dags-given-ordering" = character(0),
  "dags-with-k-roots" = "k",
  "dags-with-r-arcs" = "r"
)

#-- data, network and node types ----------------------------------------------#
available.fitted = c("bn.fit.dnet", "bn.fit.onet", "bn.fit.donet",
  "bn.fit.gnet", "bn.fit.cgnet")
available.classifiers = c("bn.naive", "bn.tan")

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

#-- graphviz plots option lists -----------------------------------------------#
graphviz.layouts = c("dot", "neato", "twopi", "circo", "fdp")
graphviz.node.shapes = c("ellipse", "circle", "rectangle")

graphviz.network.diff.methods = c("none", "from-first")

graphviz.network.diff.extra.args = list(
  "none" = character(0),
  "from-first" = c("tp.col", "tp.lty", "tp.lwd", "fp.col", "fp.lty", "fp.lwd",
                   "fn.col", "fn.lty", "fn.lwd", "show.first")
)

#-- functions to manipulate the global test counter from R --------------------#
# global test counter.
reset.test.counter = function() {

  invisible(.Call(call_reset_test_counter))

}#RESET.TEST.COUNTER

increment.test.counter = function(i = 1) {

  if (!is.real.number(i))
    stop("the increment must be a single real number.")

  invisible(.Call(call_increment_test_counter, i))

}#INCREMENT.TEST.COUNTER

test.counter = function() {

  return(.Call(call_get_test_counter))

}#TEST.COUNTER

