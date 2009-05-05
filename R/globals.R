
# Global variables.
available.discrete.tests = c("mi", "fmi", "aict", "x2", "mc-mi", "mc-x2")
available.continuous.tests = c("cor", "zf", "mi-g", "mc-mi-g", "mc-cor", "mc-zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

resampling.tests = c("mc-mi", "mc-x2", "mc-mi-g", "mc-cor", "mc-zf")

available.discrete.scores = c("lik", "loglik", "aic", "bic", "bde", "k2")
available.continuous.scores = c("bge")
available.scores = c(available.discrete.scores, available.continuous.scores)

score.equivalent.scores = c("lik", "loglik", "aic", "bic", "bde", "bge")

constraint.based.algorithms = c("gs", "iamb", "fast-iamb", "inter-iamb", "mmpc")
score.based.algorithms = c("hc")

method.labels = c(
  'gs' = "grow-shrink",
  'iamb' = "incremental association",
  'fast-iamb' = "fast incremental association",
  'inter-iamb' = "interleaved incremental association",
  'rnd' = "random/generated",
  'hc' = 'hill-climbing',
  'mmpc' = 'max-min parent children'
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
  'lik' = "likelihood",
  'loglik' = "log-likelihood",
  'bge' = "bayesian-gaussian (score equivalent)"
)

score.extra.args = list(
  "k2" = character(0),
  "bde" = "iss",
  "aic" = "k",
  "bic" = "k",
  "bge" = c("iss", "phi"),
  "lik" = character(0),
  "loglik" = character(0)
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

graphviz.enabled = FALSE
