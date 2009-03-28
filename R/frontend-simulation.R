
# generate random data from a bayesian network.
rbn = function(x, n, data, debug = FALSE) {

  # check x's class.
  check.bn(x)
  # the original data set is needed.
  check.data(data)
  # check the network against the data
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # no simulation if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # only discrete networks are supported.
  if (!is.data.discrete(data))
    stop("simulation for continuous networks not implemented.")

  rbn.discrete(x = x, n = n, data = data, debug = debug)

}#RBN

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

