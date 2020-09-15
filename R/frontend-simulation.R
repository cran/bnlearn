
# generate random data from a bayesian network.
rbn = function(x, n = 1, ...) {

  UseMethod("rbn")

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

  # call the backend.
  rbn.backend(x = x, n = n, debug = debug)

}#RBN.BN.FIT

# catch-all, fallback method.
rbn.default = function(x, n = 1, ...) {

  check.bn.or.fit(x)

}#RBN.DEFAULT

# generate a random graph.
random.graph = function(nodes, num = 1, method = "ordered", ..., debug = FALSE) {

  # check the generation method.
  check.label(method, choices = graph.generation.algorithms,
    labels = graph.generation.labels, argname = "graph generation method",
    see = "random.graph")
  # check the node labels.
  check.nodes(nodes)
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

# compute conditional probability distributions
cpdist = function(fitted, nodes, evidence, cluster = NULL, method = "ls", ...,
    debug = FALSE) {

  # check fitted's class.
  check.fit(fitted)
  # check the node labels.
  check.nodes(nodes, graph = fitted)
  # check debug.
  check.logical(debug)
  # check the generation method.
  check.label(method, choices = cpq.algorithms, labels = cpq.labels,
    argname = "sampling method", see = "cpdist")
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  extra.args = check.cpq.args(fitted = fitted, event = nodes,
                 extra.args = list(...), method = method, action = "cpdist")

  if (method == "ls") {

    if (missing(evidence))
      stop("the expression describing the evidence is missing.")

    # deparse evidence expression before passing it to the backend and beyond.
    evidence = substitute(evidence)
    # recheck event and evidence expression after deparsing.
    if (!(is.language(evidence) || identical(evidence, TRUE)))
      stop("evidence must be an unevaluated expression or TRUE.")

  }#THEN
  else if (method == "lw") {

    evidence = check.evidence(evidence, fitted)

  }#THEN

  # pass the nodes we are querying as the event of interest.
  conditional.distribution(fitted = fitted, nodes = nodes, evidence = evidence,
    method = method, extra = extra.args, cluster = cluster, debug = debug)

}#CPDIST

