
# Global variables.
available.discrete.tests = c("mh", "mi", "fmi", "aict")
available.continuous.tests = c("cor", "zf")
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
  'mh' = "Mantel-Haenszel chi-squared test",
  'mi' = "mutual information",
  'fmi' = "fast mutual information",
  'aict'= "AIC-like test",
  'cor' = "linear correlation",
  'zf' = "Fisher's Z test",
  'none' = "none"
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

# Grow-Shrink frontend.
gs = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, debug = debug,
    optimized = optimized, strict = strict, direction = direction)

}#GS

# Incremental Association frontend.
iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, method = "iamb",
    debug = debug, optimized = optimized, strict = strict,
    direction = direction)

}#IAMB

# Fast-IAMB frontend.
fast.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha,
    method = "fast-iamb", debug = debug, optimized = optimized,
    strict = strict, direction = direction)

}#FAST.IAMB

# Inter-IAMB frontend.
inter.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = NULL, alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha,
    method = "inter-iamb", debug = debug, optimized = optimized,
    strict = strict, direction = direction)

}#INTER.IAMB

# Parameter sanitization for the constraint-based learning algorithms.
# ok, it's not really a frontend (i.e. it's not exported) but it belongs here.
bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, method = "gs", debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  res = NULL
  available.methods = c("gs", "iamb", "fast-iamb", "inter-iamb")
  supported.clusters = c("MPIcluster", "PVMcluster","SOCKcluster")
  cluster.aware = FALSE

  # check the data are there.
  check.data(x)
  # check the algorithm.
  if (!(method %in% available.methods))
    stop(paste("valid values for method are:",
           paste(available.methods, collapse = " ")))
  # check test labels.
  test = check.test(test, x)
  # check the logical flags (debug, strict, optimized, direction).
  check.logical(debug)
  check.logical(strict)
  check.logical(optimized)
  check.logical(direction)
  # check alpha.
  alpha = check.alpha(alpha)

  # check cluster.
  if (!is.null(cluster)) {

    if (!(any(class(cluster) %in% supported.clusters)))
      stop("cluster is not a valid cluster object.")
    else if (!(require(snow)))
      stop("Can't find required packages: snow")
    else if (!isClusterRunning(cluster))
      stop("the cluster is stopped.")
    else {

      # enter in cluster-aware mode.
      cluster.aware = TRUE
      # set the test counter in all the cluster nodes.
      clusterEvalQ(cluster, assign(".test.counter", 0, envir = .GlobalEnv))
      # disable debugging, the slaves do not cat() here.
      if (debug) {

        warning("disabling debugging output in cluster-aware mode.")
        debug = FALSE

      }#THEN

    }#ELSE

  }#THEN

  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))

  # call the right backend.
  if (method == "gs") {

    if (cluster.aware) {

      res = grow.shrink.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = grow.shrink.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha,
        strict = strict, direction = direction, debug = debug)

    }#THEN
    else {

      res = grow.shrink(x = x, whitelist = whitelist, blacklist = blacklist,
        test = test, alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#ELSE

  }#THEN
  else if (method == "iamb") {

    if (cluster.aware) {

      res = incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#THEN
    else {

      res = incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "fast-iamb") {

    if (cluster.aware) {

      res = fast.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = fast.incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#THEN
    else {

      res = fast.incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#ELSE

  }#THEN
  else if (method == "inter-iamb") {

    if (cluster.aware) {

      res = inter.incremental.association.cluster(x = x, cluster = cluster,
        whitelist = whitelist, blacklist = blacklist, test = test,
        alpha = alpha, strict = strict, direction = direction,
        debug = debug)

    }#THEN
    else if (optimized) {

      res = inter.incremental.association.optimized(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#THEN
    else {

      res = inter.incremental.association(x = x, whitelist = whitelist,
        blacklist = blacklist, test = test, alpha = alpha, strict = strict,
        direction = direction, debug = debug)

    }#ELSE

  }#THEN

  # add tests performed by the slaves to the test counter.
  if (cluster.aware)
    res$ntests = res$ntests +
      sum(unlist(clusterEvalQ(cluster, get(".test.counter", envir = .GlobalEnv))))
  # save the learning method used.
  res$learning$algo = method
  # set the class of the return value.
  class(res) = "bn"

  invisible(res)

}#BNLEARN

# return the markov blanket of a node.
mb = function(x, node, rebuild = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # check rebuild.
  check.logical(rebuild)

  if (rebuild)
    mb.backend(x$arcs, node)
  else
    x$nodes[[node]]$mb

}#MB

# return the neighbourhood of a node.
nbr = function(x, node, rebuild = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # check rebuild.
  check.logical(rebuild)

  if (rebuild)
    nbr.backend(x$arcs, node)
  else
    x$nodes[[node]]$nbr

}#NBR

# return the arcs in the graph.
arcs = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  x$arcs

}#ARCS

# rebuild the network structure using a new set fo arcs.
"arcs<-" <- function(x, debug = FALSE, value) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a set of arcs is needed.
  if (missing(value))
    stop("no arc specified.")
  # sanitize the set of arcs.
  if (class(value) %in% c("matrix", "data.frame")) {

     if (dim(value)[2] != 2)
       stop("the arcs must have two columns.")

     if (is.data.frame(value))
       value = as.matrix(cbind(as.character(value[,1]),
         as.character(value[,2])))

  }#THEN
  else if (is.character(value)) {

    if (length(value) != 2)
      stop("the arcs must have two columns.")

    value = matrix(value, ncol = 2, byrow = TRUE,
              dimnames = list(c(), c("from", "to")))

  }#THEN
  else {

     stop("the arcs must be a matrix or data.frame with two columns.")

  }#ELSE

  # update the arcs of the network.
  x$arcs = value
  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  x

}#ARCS<-

# return the nodes in the graph.
nodes = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  names(x$nodes)

}#NODES

# build an adjacency matrix from a graph.
amat = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  arcs2amat(x$arcs, names(x$nodes))

}#AMAT

# rebuild the network structure using a new adjacency matrix.
"amat<-" <- function(x, debug = FALSE, value) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a node is needed.
  if (missing(value))
    stop("no adjacency matrix specified.")
  # check the adjacency matrix.
  check.amat(amat = value, nodes = names(x$names))

  # update the arcs of the network.
  x$arcs = amat2arcs(value, names(x$nodes))
  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  x

}#AMAT<-

# get the parents of a node.
parents = function(x, node, rebuild = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # check rebuild.
  check.logical(rebuild)

  if (rebuild)
    parents.backend(x$arcs, node)
  else
    x$nodes[[node]]$parents

}#PARENTS

# add one or more parents to a node.
"parents<-" <- function(x, node, debug = FALSE, value) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # at least one parent node is needed.
  if (missing(value))
    stop("no parent specified.")
  # node must be a valid node label.
  if (!any(value %in% names(x$nodes)))
    stop("node not present in the graph.")

  # remove duplicate labels from value.
  value = unique(value)
  # drop the parents which are not listed for inclusion.
  to.be.dropped = x$nodes[[node]]$parents[!(x$nodes[[node]]$parents %in% value)]
  # add only the nodes that were not already there.
  to.be.added = value[!(value %in% x$nodes[[node]]$parents)]

  if (debug) {

    cat("* resetting the parents of node", node, ".\n")
    cat("  > old parents: '", x$nodes[[node]]$parents, "'\n")
    cat("  > new parents: '", value, "'\n")
    cat("  > to be really dropped: '", to.be.dropped, "'\n")
    cat("  > to be really added: '", to.be.added, "'\n")

  }#THEN

  # dropping!
  for (p in to.be.dropped) {

    x = arc.operations(x = x, from = p, to = node, op = "drop",
      check.cycles = FALSE, update = FALSE, debug = debug)

  }#FOR

  # adding!
  for (p in to.be.added) {

    x = arc.operations(x = x, from = p, to = node, op = "set",
      check.cycles = TRUE, update = FALSE, debug = debug)

  }#FOR

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  x

}#PARENTS<-

# get the children of a node.
children = function(x, node, rebuild = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # check rebuild.
  check.logical(rebuild)

  if (rebuild)
    children.backend(x$arcs, node)
  else
    x$nodes[[node]]$children

}#CHILDREN

# add one or more children to a node.
"children<-" <- function(x, node, debug = FALSE, value) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = node, graph = x)
  # a node is needed.
  if (missing(value))
    stop("no children specified.")
  # node must be a valid node label.
  if (!any(value %in% names(x$nodes)))
    stop("node not present in the graph.")

  # remove duplicate labels from value.
  value = unique(value)
  # drop the parents which are not listed for inclusion.
  to.be.dropped = x$nodes[[node]]$children[!(x$nodes[[node]]$children %in% value)]
  # add only the nodes that were not already there.
  to.be.added = value[!(value %in% x$nodes[[node]]$children)]

  if (debug) {

    cat("* resetting the children of node", node, ".\n")
    cat("  > old children: '", x$nodes[[node]]$children, "'\n")
    cat("  > new children: '", value, "'\n")
    cat("  > to be really dropped: '", to.be.dropped, "'\n")
    cat("  > to be really added: '", to.be.added, "'\n")

  }#THEN

  # dropping!
  for (child in to.be.dropped) {

    x = arc.operations(x = x, from = node, to = child, op = "drop",
      check.cycles = FALSE, update = FALSE, debug = debug)

  }#FOR

  # adding!
  for (child in to.be.added) {

    x = arc.operations(x = x, from = node, to = child, op = "set",
      check.cycles = TRUE, update = FALSE, debug = debug)

  }#FOR

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  x

}#CHILDREN<-

# get the root nodes of a bayesian network.
rootnodes = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  rootnodes.backend(x$arcs, names(x$nodes))

}#ROOTNODES

# get the leaf nodes of a bayesian network.
leafnodes = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  leafnodes.backend(x$arcs, names(x$nodes))

}#LEAFNODES

# get the number of paraters of the bayesian network.
nparams = function(x, data, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check the data are there.
  check.data(data)
  # only discrete bayesian networks are supported.
  if (!is.data.discrete(data))
    stop("parameter enumeration for continuous networks not implemented.")
  # check debug.
  check.logical(debug)
  # nparams is unknown for partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  params = nparams.backend(x, data, real = TRUE)

  if (debug) {

    cat(">  number of parameters per node:\n")
    print(params)

  }#THEN

  sum(params)

}#NPARAMS

# check if a graph is acyclic.
acyclic = function(x, directed, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check debug.
  check.logical(debug)

  if (missing(directed)) {

    is.acyclic(x$arcs, names(x$nodes), debug = debug)

  }#THEN
  else {

    # check directed.
    check.logical(directed)

    if (directed)
      is.dag.acyclic(arcs = x$arcs, nodes = names(x$nodes), debug = debug)
    else
      is.pdag.acyclic(arcs = x$arcs, nodes = names(x$nodes), debug = debug)

  }#ELSE

}#ACYCLIC

# check if there's a path between two specific nodes.
path = function(x, from, to, direct = TRUE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a valid node is needed.
  check.node(node = from, graph = x)
  # another valid node is needed.
  check.node(node = to, graph = x)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")

  if (is.pdag(x$arcs, names(x$nodes)))
    has.pdag.path(from, to, names(x$nodes),
      arcs2amat(x$arcs, names(x$nodes)),
      exclude.direct = !direct)
  else
    has.path(from, to, names(x$nodes),
      arcs2amat(x$arcs, names(x$nodes)),
      exclude.direct = !direct)

}#PATH

# describe the network with a "model string".
modelstring = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # no model string if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  formula.backend(x)

}#MODELSTRING

# return the partial node ordering implied by the graph structure.
node.ordering = function(x, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check debug.
  check.logical(debug)
  # no model string if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  schedule(x, debug = debug)

}#NODE.ORDERING

# generate an object of class bn from a model string.
model2network = function(string, debug = FALSE) {

  # check string's class.
  if (!is(string, "character"))
    stop("string must be a character string.")
  # check debug.
  check.logical(debug)

  model2network.backend(string, debug = debug)

}#MODEL2NETWORK

# generate a valid blacklist from a partial node ordering.
ordering2blacklist = function(nodes) {

  # nodes must be a vector of character strings.
  if (!is(nodes, "character"))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # at least one node is needed.
  if (length(nodes) < 1)
    stop("at leat one node label is needed.")
  # no duplicates allowed.
  if(any(duplicated(nodes)))
     stop("node labels must be unique.")

  do.call(rbind,

    sapply(seq(from = 1, to = length(nodes)),
      function(i) {

        cbind(from = rep(nodes[i], i - 1), to = nodes[(1:i) - 1])

    })
  )

}#ORDERING2BLACKLIST

# generate random data from a network.
rbn = function(x, n, data, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # the original data set is needed.
  check.data(data)
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

# set an arc direction manually.
set.arc = function(x, from, to, check.cycles = TRUE, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "set",
    check.cycles = check.cycles, debug = debug)

}#SET.ARC

# remove an arc from the graph.
drop.arc = function(x, from, to, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "drop",
    check.cycles = FALSE, debug = debug)

}#DROP.ARC

# reverse an arc in the graph.
reverse.arc = function(x, from, to, check.cycles = TRUE, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "reverse",
    check.cycles = check.cycles, debug = debug)

}#REVERSE.ARC

# Parameter sanitization for arc operations.
# ok, it's not really a frontend (i.e. it's not exported) but it belongs here.
arc.operations = function(x, from, to, op = NULL, check.cycles, update = TRUE, debug = FALSE) {

available.ops = c("set", "drop", "reverse")

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check the op code.
  if (!(op %in% available.ops))
    stop("valid op codes are 'set', 'drop' and 'reverse'.")
  # a valid node is needed.
  check.node(node = from, graph = x)
  # another valid node is needed.
  check.node(node = to, graph = x)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check logical flags (debug, check.cycles, update).
  check.logical(debug)
  check.logical(check.cycles)
  check.logical(update)

  # add/reverse/orient the arc.
  if (op == "set") {

    if (debug) cat("* setting arc", from, "->", to, ".\n")
    x$arcs = set.arc.direction(from, to, x$arcs, debug = debug)

  }#THEN
  else if (op == "drop") {

    if (debug) cat("* dropping any arc between ", from, "and", to, ".\n")
    x$arcs = drop.arc.backend(x$arcs, c(from, to), debug = debug)

  }#THEN
  else if (op == "reverse") {

    if (debug) cat("* reversing any arc between ", from, "and", to, ".\n")
    x$arcs = reverse.arc.backend(from, to, x$arcs, debug = debug)

  }#THEN

  # check whether the graph is still acyclic; not needed if an arc is dropped.
  if (check.cycles && (op != "drop"))
    if (!is.acyclic(x$arcs, names(x$nodes)))
      stop("the resulting graph contains cycles.")

  # update the network structure.
  if (update)
    x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  invisible(x)

}#ARC.OPERATIONS

# compute the score of a network.
score = function(x, data, type = NULL, ..., debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # the original data set is needed.
  check.data(data)
  # check debug.
  check.logical(debug)
  # no score if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the score label.
  type = check.score(type, data)

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = type, network = x, 
                 data = data, extra.args = list(...))

  # compute the network score.
  network.score(network = x, data = data, score = type,
    extra.args = extra.args, debug = debug)

}#SCORE

# create an empty graph from a given set of nodes.
empty.graph = function(nodes) {

  # nodes must be a vector of character strings.
  if (!is(nodes, "character"))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # at least one node is needed.
  if (length(nodes) < 1)
    stop("at leat one node label is needed.")
  # no duplicates allowed.
  if(any(duplicated(nodes)))
     stop("node labels must be unique.")

  empty.graph.backend(nodes)

}#EMPTY.GRAPH

# generate a random graph.
random.graph = function(nodes, prob = 0.5) {

  # nodes must be a vector of character strings.
  if (!is(nodes, "character"))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # at least one node is needed.
  if (length(nodes) < 1)
    stop("at leat one node label is needed.")
  # no duplicates allowed.
  if (any(duplicated(nodes)))
    stop("node labels must be unique.")
  # prob must be numeric.
  if (!is(prob, "numeric"))
    stop("the branching probability must be a numeric value.")
  # prob must be scalar.
  if (length(prob) > 1)
    stop("the branching probability must be a scalar.")
  # prob must be a value in [0,1].
  if ((prob < 0) || (prob > 1))
    stop("the branching probability must be a numeric value in [0,1].")

  random.graph.backend(nodes = nodes, prob = prob)

}#RANDOM.GRAPH

# compare two bayesian network structures.
compare = function (r1, r2, debug = FALSE) {

  result = TRUE

  # check both objects' class.
  if (!is(r1, "bn") || !is(r2, "bn"))
    stop("both r1 and r2 must be objects of class 'bn'.")
  # check debug.
  check.logical(debug)

  # check the two graphs have the same nodes.
  r1.nodes = names(r1$nodes)
  r2.nodes = names(r2$nodes)

  if (!identical(sort(r1.nodes), sort(r2.nodes))) {

    if (debug) {

      cat("* nodes in r1 not present in r2:\n")
      print(r1.nodes[!(r1.nodes %in% r2.nodes)])
      cat("* nodes in r2 not present in r1:\n")
      print(r2.nodes[!(r2.nodes %in% r1.nodes)])

    }#THEN

    return(FALSE)

  }#THEN

  # for each node check ...
  check = sapply(names(r1$nodes),

    function(node) {

      node.result = TRUE
      r1.node = r1$nodes[[node]]
      r2.node = r2$nodes[[node]]

      # ... the markov blanket ...
      if (!identical(sort(r1.node$mb), sort(r2.node$mb))) {

        if (debug) {

          cat("* nodes in the markov blanket of", node, "in r1 not present in r2:\n")
          print(r1.node$mb[!(r1.node$mb %in% r2.node$mb)])
          cat("* nodes in the markov blanket of", node, "in r2 not present in r1:\n")
          print(r2.node$mb[!(r2.node$mb %in% r1.node$mb)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... and the neighbourhood ...
      if (!identical(sort(r1.node$nbr), sort(r2.node$nbr))) {

        if (debug) {

          cat("* nodes in the neighbourhood of", node, "in r1 not present in r2:\n")
          print(r1.node$nbr[!(r1.node$nbr %in% r2.node$nbr)])
          cat("* nodes in the neighbourhood of", node, "in r2 not present in r1:\n")
          print(r2.node$nbr[!(r2.node$nbr %in% r1.node$nbr)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... the parents ...
      if (!identical(sort(r1.node$parents), sort(r2.node$parents))) {

        if (debug) {

          cat("* parents of", node, "in r1 not present in r2:\n")
          print(r1.node$parents[!(r1.node$parents %in% r2.node$parents)])
          cat("* parents of", node, "in r2 not present in r1:\n")
          print(r2.node$parents[!(r2.node$parents %in% r1.node$parents)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... and the children.
      if (!identical(sort(r1.node$children), sort(r2.node$children))) {

        if (debug) {

          cat("* children of", node, "in r1 not present in r2:\n")
          print(r1.node$children[!(r1.node$children %in% r2.node$children)])
          cat("* children of", node, "in r2 not present in r1:\n")
          print(r2.node$children[!(r2.node$children %in% r1.node$children)])

        }#THEN

        node.result = FALSE

      }#THEN

      return(node.result)

    }

  )

  if (!all(check)) result = FALSE

  # check directed arcs.
  r1.arcs = apply(r1$arcs[!is.undirected(r1$arcs), , drop = FALSE], 1, paste, collapse = " -> ")
  r2.arcs = apply(r2$arcs[!is.undirected(r2$arcs), , drop = FALSE], 1, paste, collapse = " -> ")

  if (!identical(sort(r1.arcs), sort(r2.arcs))) {

    if (debug) {

      cat("* directed arcs in r1 not present in r2:\n")
      print(r1.arcs[!(r1.arcs %in% r2.arcs)])
      cat("* directed arcs in r2 not present in r1:\n")
      print(r2.arcs[!(r2.arcs %in% r1.arcs)])

    }#THEN

    result = FALSE

  }#THEN

  # check undirected arcs.
  r1.arcs = apply(r1$arcs[is.undirected(r1$arcs), , drop = FALSE], 1, paste, collapse = " - ")
  r2.arcs = apply(r2$arcs[is.undirected(r2$arcs), , drop = FALSE], 1, paste, collapse = " - ")

  if (!identical(sort(r1.arcs), sort(r2.arcs))) {

    if (debug) {

      cat("* undirected arcs in r1 not present in r2:\n")
      print(r1.arcs[!(r1.arcs %in% r2.arcs)])
      cat("* undirected arcs in r2 not present in r1:\n")
      print(r2.arcs[!(r2.arcs %in% r1.arcs)])

    }#THEN

    result = FALSE

  }#THEN

  result

}#COMPARE

# infer the direction of an ipothetic arc between two specified nodes.
choose.direction = function(x, arc, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check the data are there.
  check.data(data)
  # check the arc is there.
  check.arc(arc, x)
  # check debug.
  check.logical(debug)
  # check criterion.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data)
    else
      criterion = x$learning$test

  }#THEN
  else {

    if (criterion %in% available.tests)
      criterion = check.test(criterion, data)
    else if (criterion %in% available.scores)
      criterion = check.score(criterion, data)
    else stop(paste("valid criteria are:",
           paste(available.tests, available.scores, collapse = " ")))

  }#ELSE

  if (debug)
    cat("* testing", arc[1], "-", arc[2], "for direction.\n" )

  if (criterion %in% available.tests) {

    alpha = check.alpha(list(...)$alpha, network = x)

    x = choose.direction.test(x, data = data, arc = arc, test = criterion,
          alpha = alpha, debug = debug)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x, 
                   data = data, extra.args = list(...))

    x = choose.direction.score(x, data = data, arc = arc, score = criterion,
          extra.args = extra.args, debug = debug)

  }#ELSE

  invisible(x)

}#CHOOSE.DIRECTION

# Hill Climbing greedy search frontend.
hc = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = NULL, ..., debug = FALSE, restart = 0, perturb = 1,
    optimized = TRUE) {

  greedy.search(x = x, start = start, whitelist = whitelist,
    blacklist = blacklist, score = score, heuristic = "hc",
    ..., debug = debug, restart = restart, perturb = perturb,
   optimized = optimized)

}#HC

# Parameter sanitization for the score-based learning algorithms.
# ok, it's not really a frontend (i.e. it's not exported) but it belongs here.
greedy.search = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
    score = "k2", heuristic = "hc", ..., debug = FALSE, restart = 0,
    perturb = 1, optimized = FALSE) {

  # check the data are there.
  check.data(x)
  # check the score label.
  score = check.score(score, x)
  # check debug.
  check.logical(debug)
  # check restart and perturb.
  check.restart(restart, perturb)
  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = names(x))
  else {

    # check start's class.
    if (!is(start, "bn"))
      stop("x must be an object of class 'bn'.")
    # set all nodes as updated if the preseed network is not empty, 
    # so that all cache lookups are skipped.
    if (nrow(start$arcs) > 0)
      start$updates = array(rep(0, length(start$nodes)), 
                        dimnames = list(names(start$nodes)))

  }#ELSE

  # apply the whitelist to the preseeded network.
  if (!is.null(whitelist)) {

    for (i in 1:nrow(whitelist))
      start$arcs = set.arc.direction(whitelist[i, "from"],
                       whitelist[i, "to"], start$arcs)

  }#THEN

  # apply the blacklist to the preseeded network.
  if (!is.null(blacklist)) {

    blacklisted = apply(start$arcs, 1, function(x){ is.blacklisted(blacklist, x) })
    start$arcs = start$arcs[!blacklisted, , drop = FALSE]

  }#THEN

  # be sure the graph structure is up to date.
  start$nodes = cache.structure(names(start$nodes), start$arcs)
  # no party if the graph is partially directed.
  if (is.pdag(start$arcs, names(start$nodes)))
    stop("the graph is only partially directed.")
  # check whether the graph is acyclic.
  if (!is.dag.acyclic(start$arcs, names(start$nodes)))
    stop("the preseeded graph contains cycles.")

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = start, 
                 data = x, extra.args = list(...))

  # create the test counter in .GlobalEnv.
  assign(".test.counter", 0, envir = .GlobalEnv)

  if (heuristic == "hc") {

    if (optimized) {

      res = hill.climbing.optimized(x = x, start = start, whitelist = whitelist,
        blacklist = blacklist, score = score, extra.args = extra.args,
        restart = restart, perturb = perturb, debug = debug)

    }#THEN
    else {

      res = hill.climbing(x = x, start = start, whitelist = whitelist,
        blacklist = blacklist, score = score, extra.args = extra.args,
        restart = restart, perturb = perturb, debug = debug)

    }#ELSE

  }#THEN

  # set the metadata of the network.
  res$learning$algo = heuristic
  res$learning$ntests = get(".test.counter", envir = .GlobalEnv)
  res$learning$test = score
  res$learning$args = extra.args

  invisible(res)

}#GREEDY.SEARCH

# measure the strength of the arcs in a directed graph.
arc.strength = function(x, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the data are there.
  check.data(data)
  # check debug.
  check.logical(debug)
  # check criterion.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data)
    else
      criterion = x$learning$test

  }#THEN
  else {

    if (criterion %in% available.tests)
      criterion = check.test(criterion, data)
    else if (criterion %in% available.scores)
      criterion = check.score(criterion, data)
    else stop(paste("valid criteria are:",
           paste(available.tests, available.scores, collapse = " ")))

  }#ELSE

  # expand and sanitize score-specific arguments and the alpha threshold.
  if (criterion %in% available.tests) {

    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize the alpha threshold.
    arc.strength.test(network = x, data = data, alpha = alpha, 
      test = criterion, debug = debug)

  }#THEN
  else {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x, 
                 data = data, extra.args = list(...))

    arc.strength.score(network = x, data = data, score = criterion, 
      extra = extra.args, debug = debug)

  }#ELSE

}#ARC.STRENGTH

