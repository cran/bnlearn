
# generate random data from a bayesian network.
rbn = function(x, n = 1, ...) {

  UseMethod("rbn", x)

}#RBN

# generate random data from a network structure (parameters are learned on the fly).
rbn.bn = function(x, n = 1, data, fit = "mle", ..., debug = FALSE) {

  # fit the parameters of the bayesian network.
  fitted = bn.fit(x = x, data = data, method = fit, ..., debug = debug)
  # call the other method.
  rbn.bn.fit(x = fitted, n = n, debug = debug)

}#RBN.BN

# generate random data from a bayesian network.
rbn.bn.fit = function(x, n = 1, ..., debug = FALSE) {

  # check the size of the sample to be generated.
  if (!is.positive.integer(n))
    stop("the number of observations to be generated must be a positive integer number.")
  # check debug.
  check.logical(debug)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # call the right backend.
  if (is.fitted.discrete(x))
    rbn.discrete(x = x, n = n, data = NULL, debug = debug)
  else
    rbn.continuous(x = x, n = n, data = NULL, debug = debug)

}#RBN.BN.FIT

# catch-all, fallback method.
rbn.default = function(x, n = 1, ...) {

  check.bn.or.fit(x)

}#RBN.DEFAULT

# generate a random graph.
random.graph = function(nodes, num = 1, method = "ordered", ..., debug = FALSE) {

  # check the generation method.
  if (!(method %in% graph.generation.algorithms))
    stop(paste("valid generation methods are:",
           paste(graph.generation.algorithms, collapse = " ")))
  # check the node labels.
  check.nodes(nodes, min.nodes = 3)
  # check the number of graph to generate.
  if (!is.positive.integer(num))
    stop(" the number of graphs to generate must be a positive integer number.")

  # expand and sanitize method-specific arguments.
  extra.args = check.graph.generation.args(method = method,
                 nodes = nodes, extra.args = list(...))

  random.graph.backend(num = num, nodes = nodes, method = method,
    extra.args = extra.args, debug = debug)

}#RANDOM.GRAPH

# create an empty graph from a given set of nodes.
empty.graph = function(nodes, num = 1) {

  random.graph(nodes = nodes, num = num, method = "empty", debug = FALSE)

}#EMPTY.GRAPH

# perform conditional probability queries.
cpquery = function(fitted, event, evidence, cluster = NULL, method = "ls", ..., debug = FALSE) {

  # check fitted's class.
  check.fit(fitted)
  # check debug.
  check.logical(debug)
  # check event and evidence.
  if (missing(event))
    stop("the expression describing the event is missing.")
  if (missing(evidence))
    stop("the expression describing the evidence is missing.")
  # check the generation method.
  if (!(method %in% cpq.algorithms))
    stop(paste(c("valid conditional probability query methods are:\n",
      sprintf("    %-10s %s\n", cpq.algorithms, cpq.labels[cpq.algorithms])), sep = ""))
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output in cluster-aware mode.")
      debug = FALSE

    }#THEN

  }#THEN

  extra.args = check.cpq.args(fitted = fitted, extra.args = list(...),
                 method = method)

  # deparse both event and evidence expressions before passing them to
  # the backend and beyond.
  event = substitute(event)
  evidence = substitute(evidence)

  # recheck event and evidence expression after deparsing.
  if (!(is.language(event) || identical(event, TRUE)))
    stop("event must be an unevaluated expression or TRUE.")
  if (!(is.language(evidence) || identical(evidence, TRUE)))
    stop("evidence must be an unevaluated expression or TRUE.")

  conditional.probability.query(fitted = fitted, event = event,
    evidence = evidence, method = method, extra = extra.args,
    probability = TRUE, cluster = cluster, debug = debug)

}#CPQUERY

# compute conditional probability distributions
cpdist = function(fitted, nodes, evidence, cluster = NULL, method = "ls", ..., debug = FALSE) {

  # check fitted's class.
  check.fit(fitted)
  # check the node labels.
  check.nodes(nodes, graph = fitted, min.nodes = 1)
  # check debug.
  check.logical(debug)
  if (missing(evidence))
    stop("the expression describing the evidence is missing.")
  # check the generation method.
  if (!(method %in% cpq.algorithms))
    stop(paste(c("valid conditional probability query methods are:\n",
      sprintf("    %-10s %s\n", cpq.algorithms, cpq.labels[cpq.algorithms])), sep = ""))
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output in cluster-aware mode.")
      debug = FALSE

    }#THEN

  }#THEN

  extra.args = check.cpq.args(fitted = fitted, extra.args = list(...),
                 method = method)

  # deparse evidence expression before passing it to the backend and beyond.
  evidence = substitute(evidence)
  # recheck event and evidence expression after deparsing.
  if (!(is.language(evidence) || identical(evidence, TRUE)))
    stop("evidence must be an unevaluated expression or TRUE.")

  # pass the nodes we are querying as the event of interest.
  conditional.probability.query(fitted = fitted, event = nodes,
    evidence = evidence, method = method, extra = extra.args,
    probability = FALSE, cluster = cluster, debug = debug)

}#CPDIST

