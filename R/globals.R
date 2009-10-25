
# Global variables.
available.discrete.tests = c("mi", "fmi", "aict", "x2", "mc-mi", "mc-x2")
available.continuous.tests = c("cor", "zf", "mi-g", "mc-mi-g", "mc-cor", "mc-zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

resampling.tests = c("mc-mi", "mc-x2", "mc-mi-g", "mc-cor", "mc-zf")

available.discrete.scores = c("loglik", "aic", "bic", "bde", "k2")
available.continuous.scores = c("bge", "loglik-g", "aic-g", "bic-g")
available.scores = c(available.discrete.scores, available.continuous.scores)

score.equivalent.scores = c("loglik", "aic", "bic", "bde", "bge", "loglik-g", "aic-g", "bic-g")

constraint.based.algorithms = c("gs", "iamb", "fast.iamb", "inter.iamb", "mmpc")
score.based.algorithms = c("hc")
hybrid.algorithms = c("rshc", "mmhc")
available.learning.algorithms = c(constraint.based.algorithms, score.based.algorithms, hybrid.algorithms)

always.dag.result = c(score.based.algorithms, hybrid.algorithms)

available.mvber.vartests = c("tvar", "gvar", "nvar")

method.labels = c(
  'gs' = "grow-shrink",
  'iamb' = "incremental association",
  'fast.iamb' = "fast incremental association",
  'inter.iamb' = "interleaved incremental association",
  'rnd' = "random/generated",
  'hc' = 'hill-climbing',
  'mmpc' = 'max-min parent children',
  'rshc' = 'restricted hill climbing',
  'mmhc' = 'max-min hill climbing'
)

test.labels = c(
  'mi' = "mutual information (discrete)",
  'mc-mi' = "mutual information (discrete, Monte Carlo)",
  'mi-g' = "mutual information (gaussian)",
  'mc-mi-g' = "mutual information (gaussian, Monte Carlo)",
  'fmi' = "fast mutual information",
  'aict'= "AIC-like test",
  'x2'= "Pearson's X^2",
  'mc-x2'= "Pearson's X^2 (Monte Carlo)",
  'cor' = "linear correlation",
  'mc-cor' = "linear correlation (Monte Carlo)",
  'zf' = "Fisher's Z test",
  'mc-zf' = "Fisher's Z test (Monte Carlo)"
)

score.labels = c(
  'k2' = "Cooper & Herskovits' K2",
  'bde' = "bayesian-dirichlet (score equivalent)",
  'aic' = "Akaike information criterion",
  'bic' = "bayesian information criterion",
  'loglik' = "log-likelihood",
  'bge' = "bayesian-gaussian (score equivalent)",
  'loglik-g' = "log-likelihood (gaussian)",
  'aic-g' = "Akaike information criterion (gaussian)",
  'bic-g' = "bayesian information criterion (gaussian)"
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

graph.generation.algorithms = c("ordered", "ic-dag", "empty")

graph.generation.labels = c(
  "ordered" = "full ordering",
  "ic-dag" = "multiconnected dags",
  "empty" = "empty"
)

graph.generation.extra.args = list(
  "ordered" = "prob",
  "ic-dag" = c("burn.in", "max.degree", "max.in.degree", "max.out.degree")
)

mvber.labels = list(
  "tvar" = "total variance",
  "gvar" = "generalized variance",
  "nvar" = "squared Frobenius norm"
)

graphviz.enabled = FALSE
lattice.enabled = FALSE
