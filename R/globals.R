
# Global variables.
available.discrete.tests = c("mi", "fmi", "aict", "x2")
available.continuous.tests = c("cor", "zf", "mi-g")
available.tests = c(available.discrete.tests, available.continuous.tests)

available.discrete.scores = c("lik", "loglik", "aic", "bic", "dir", "bde", "k2")
available.continuous.scores = c("bge")
available.scores = c(available.discrete.scores, available.continuous.scores)

score.equivalent.scores = c("lik", "loglik", "aic", "bic", "dir", "bde", "bge")

method.labels = c(
  'gs' = "grow-shrink",
  'iamb' = "incremental association",
  'fast-iamb' = "fast incremental association",
  'inter-iamb' = "interleaved incremental association",
  'rnd' = "random/generated",
  'hc' = 'hill-climbing'
)

test.labels = c(
  'mi' = "mutual information (discrete)",
  'mi-g' = "mutual information (gaussian)",
  'fmi' = "fast mutual information",
  'aict'= "AIC-like test",
  'x2'= "Pearson's X^2",
  'cor' = "linear correlation",
  'zf' = "Fisher's Z test"
)

score.labels = c(
  'k2' = "Cooper & Herskovits' K2",
  'bde' = "bayesian-dirichlet (score equivalent)",
  'dir' = "bayesian-dirichlet (score equivalent)",
  'aic' = "Akaike information criterion",
  'bic' = "bayesian information criterion",
  'lik' = "likelihood",
  'loglik' = "log-likelihood",
  'bge' = "bayesian-gaussian (score equivalent)"
)

score.extra.args = list(
  "k2" = character(0),
  "dir" = "iss",
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
