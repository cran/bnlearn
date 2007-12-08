
# Global variables.
available.discrete.tests = c("mh", "mi", "fmi")
available.continuous.tests = c("cor", "zf")
available.tests = c(available.discrete.tests, available.continuous.tests)

method.labels = c(
  'gs' = "grow-shrink", 
  'iamb' = "incremental association", 
  'fast-iamb' = "fast incremental association",
  'inter-iamb' = "interleaved incremental association",
  'rnd' = "random/generated"
)

test.labels = c(
  'mh' = "Mantel-Haenszel chi-squared test",
  'mi' = "mutual information",
  'fmi' = "fast mutual information",
  'cor' = "linear correlation",
  'zf' = "Fisher's Z test"
)

# Grow-Shrink frontend.
gs = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, debug = debug,
    optimized = optimized, strict = strict, direction = direction)

}#GS

# Incremental Association frontend.
iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha, method = "iamb",
    debug = debug, optimized = optimized, strict = strict,
    direction = direction)

}#IAMB

# Fast-IAMB frontend.
fast.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha,
    method = "fast-iamb", debug = debug, optimized = optimized,
    strict = strict, direction = direction)

}#FAST.IAMB

# Inter-IAMB frontend.
inter.iamb = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mi", alpha = 0.05, debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  bnlearn(x = x, cluster = cluster, whitelist = whitelist,
    blacklist = blacklist, test = test, alpha = alpha,
    method = "inter-iamb", debug = debug, optimized = optimized,
    strict = strict, direction = direction)

}#INTER.IAMB

# Parameter sanitization for the learning algorithms' frontends.
# ok, it's not really a frontend (i.e. it's not exported) but it belongs here.
bnlearn = function(x, cluster = NULL, whitelist = NULL, blacklist = NULL,
    test = "mh", alpha = 0.05, method = "gs", debug = FALSE, optimized = TRUE,
    strict = TRUE, direction = FALSE) {

  assign(".test.counter", 0, envir = .GlobalEnv)

  res = NULL
  available.methods = c("gs", "iamb", "fast-iamb", "inter-iamb")
  supported.clusters = c("MPIcluster", "PVMcluster","SOCKcluster")
  cluster.aware = FALSE

  # check the data are there.
  if (missing(x))
    stop("the data are missing.")
  # x must be a data frame.
  if(!is.data.frame(x))
    stop("x must be a data frame.")
  # check the data for NULL/NaN/NA.
  if (missing.data(x))
    stop("the data set contains NULL/NaN/NA values.")
  # check the test.
  if (!(test %in% available.tests))
    stop(paste("valid values for test are:",
           paste(available.tests, collapse = " ")))
  # check the algorithm.
  if (!(method %in% available.methods))
    stop(paste("valid values for method are:",
           paste(available.methods, collapse = " ")))
  # check if it's the right test for the data (discrete, continuous).
  if (!is.data.discrete(x) && (test %in% available.discrete.tests))
    stop(paste("test '", test, "' may be used with discrete data only.", sep = ""))
  if (!is.data.continuous(x) && (test %in% available.continuous.tests))
    stop(paste("test '", test, "' may be used with continuous data only.", sep = ""))
  # check the variables are either all continuous or all discrete.
  if (!is.data.discrete(x) && !is.data.continuous(x))
    stop("variables must be either all numeric or all factors.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # check strict.
  if (!is.logical(strict) || is.na(strict))
    stop("strict must be a logical value (TRUE/FALSE).")
  # check optimized.
  if (!is.logical(optimized) || is.na(optimized))
    stop("optimized must be a logical value (TRUE/FALSE).")
  # check direction.
  if (!is.logical(direction) || is.na(direction))
    stop("direction must be a logical value (TRUE/FALSE).")
  # check alpha.
  if (!is.numeric(alpha) || (alpha > 1) || (alpha < 0))
    stop("alpha must be a numerical value in [0,1].")

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

  # NOTE: whitelist and blacklist relationship is the same as hosts.allow
  # and hosts.deny.

  if (!is.null(whitelist)) {

    if (class(whitelist) %in% c("matrix", "data.frame")) {

       if (dim(whitelist)[2] != 2)
         stop("whitelist must have two columns.")

       if (is.data.frame(whitelist))
         whitelist = as.matrix(cbind(as.character(whitelist[,1]),
           as.character(whitelist[,2])))

    }#THEN
    else if (is.character(whitelist)) {

      if (length(whitelist) != 2)
        stop("whitelist must have two columns.")

      whitelist = matrix(whitelist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("whitelist must be a matrix or data.frame with two columns.")

    }#ELSE

    # add column names for easy reference.
    colnames(whitelist) = c("from", "to")
    # drop duplicate rows.
    whitelist = unique(whitelist)

    # check all the names in the whitelist against the column names of x.
    if (any(!(unique(as.vector(whitelist)) %in% names(x))))
     stop("unknown node label present in the whitelist.")

    # if the whitelist itself contains cycles, no acyclic graph
    # can be learned.
    if (!is.acyclic(whitelist, colnames(x)))
      stop("this whitelist does not allow an acyclic graph.")


  }#THEN

  if (!is.null(blacklist)) {

    if (class(blacklist) %in% c("matrix", "data.frame")) {

       if (dim(blacklist)[2] != 2)
         stop("blacklist must have two columns.")

       if (is.data.frame(blacklist))
         blacklist = as.matrix(cbind(as.character(blacklist[,1]),
           as.character(blacklist[,2])))

    }#THEN
    else if (is.character(blacklist)) {

      if (length(blacklist) != 2)
        stop("blacklist must have two columns.")

      blacklist = matrix(blacklist, ncol = 2, byrow = TRUE)

    }#THEN
    else {

      stop("blacklist must be a matrix or data.frame with two columns.")

    }#ELSE

    # add column names for easy reference.
    colnames(blacklist) = c("from", "to")
    # drop duplicate rows.
    blacklist = unique(blacklist)

    # check all the names in the blacklist against the column names of x.
    if (any(!(unique(as.vector(blacklist)) %in% names(x))))
      stop("unknown node label present in the blacklist.")

  }#THEN

  # if x -> y is whitelisted but y -> x is not, it is to be blacklisted.
  if (!is.null(whitelist)) {

    apply(whitelist, 1,
      function(x) {
        if (!is.whitelisted(whitelist, x[c(2,1)]))
          assign("blacklist", rbind(blacklist, x[c(2,1)]),
            envir = sys.frame(-2))
      }#FUNCTION
    )

  }#THEN

  # if x -> y is whitelisted, it is to be removed from the blacklist.
  if (!is.null(blacklist) && !is.null(whitelist)) {

    blacklist = blacklist[!apply(blacklist, 1,
      function(x){ is.whitelisted(whitelist, x) }),]

    blacklist = matrix(blacklist, ncol = 2, byrow = FALSE,
      dimnames = list(NULL, c("from", "to")))

  }#THEN

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
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
  # check rebuild.
  if (!is.logical(rebuild) || is.na(rebuild))
    stop("rebuild must be a logical value (TRUE/FALSE).")

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
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
  # check rebuild.
  if (!is.logical(rebuild) || is.na(rebuild))
    stop("rebuild must be a logical value (TRUE/FALSE).")

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
  # the adjacency matrix must, well, be a matrix.
  if (!is(value, "matrix"))
    stop("this is not a valid adjacency matrix.")
  # column names must be a valid node label.
  if (!is.null(colnames(value)))
    if (!all(colnames(value) %in% names(x$nodes)))
      stop("node (column label) not present in the graph.")
  # column names must be a valid node label.
  if (!is.null(rownames(value)))
    if (!all(rownames(value) %in% names(x$nodes)))
      stop("node (row label) not present in the graph.")

  # update the arcs of the network.
  x$arcs = amat2arcs(value, names(x$nodes), debug = debug)
  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  x

}#AMAT<-

# get the parents of a node.
parents = function(x, node, rebuild = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
  # check rebuild.
  if (!is.logical(rebuild) || is.na(rebuild))
    stop("rebuild must be a logical value (TRUE/FALSE).")

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
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
  # a node is needed. 
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
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
  # check rebuild.
  if (!is.logical(rebuild) || is.na(rebuild))
    stop("rebuild must be a logical value (TRUE/FALSE).")

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
  # a node is needed. 
  if (missing(node))
    stop("no node specified.")
  # a node label must be a character string.
  if (!is(node, "character"))
    stop("node must be a character string, the label of a node.")
  # only one node is needed.
  if (length(node) > 1)
    stop("only a single node may be specified.")
  # node must be a valid node label.
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")
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
  # only discrete bayesian networks are supported.
  if (!is.data.discrete(data))
    stop("parameter enumeration for continuous networks not implemented.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # nparams is unknown for partially directed graphs.
  if (any(is.undirected(x$arcs)))
    stop("the graph is only partially directed.")

  params = nparams.backend(x, data, real = TRUE)

  if (debug) {

    cat(">  number of parameters per node:\n")
    print(params)

  }#THEN

  sum(params)

}#NPARAMS

acyclic = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  is.acyclic(x$arcs, names(x$nodes))

}#ACYCLIC

path = function(x, from, to, direct = TRUE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # both 'from' and 'to' are needed.
  if (missing(from) || missing(to))
    stop("'from' and/or 'to' are missing.")
  # 'from' must be a character string.
  if (!is(from, "character"))
    stop("'from' must be a character string, the label of a node.")
  # only one 'from' node is needed.
  if (length(from) > 1)
    stop("only a single 'from' node may be specified.")
  # 'to' must be a character string.
  if (!is(to, "character"))
    stop("'to' must be a character string, the label of a node.")
  # only one 'to' node is needed.
  if (length(to) > 1)
    stop("only a single 'to' node may be specified.")
  # 'from' must be a valid node label.
  if (!(from %in% names(x$nodes)))
    stop("'from' not present in the graph.")
  # 'to' must be a valid node label.
  if (!(to %in% names(x$nodes)))
    stop("'to' not present in the graph.")
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")

  has.path(from, to, names(x$nodes), 
    arcs2amat(x$arcs, names(x$nodes)), 
    exclude.direct = FALSE)

}#PATH

# describe the network with a "model string".
modelstring = function(x) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # no model string if the graph is partially directed.
  if (any(is.undirected(x$arcs)))
    stop("the graph is only partially directed.")

  formula.backend(x)

}#MODELSTRING

# return the partial node ordering implied by the graph structure.
node.ordering = function(x, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # no model string if the graph is partially directed.
  if (any(is.undirected(x$arcs)))
    stop("the graph is only partially directed.")

  schedule(x, debug = debug)

}#NODE.ORDERING

# generate an object of class bn from a model string.
model2network = function(string, debug = FALSE) {

  # check string's class.
  if (!is(string, "character"))
    stop("string must be a character string.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")

  model2network.backend(string, debug = debug)

}#MODEL2NETWORK

# generate random data from a network.
rbn = function(x, n, data, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # the original data set is needed.
  if (missing(data))
    stop("the data are missing.")
  # check the data for NULL/NaN/NA.
  if (missing.data(data))
    stop("the data set contains NULL/NaN/NA values.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # no simulation if the graph is partially directed.
  if (any(is.undirected(x$arcs)))
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

# remove an arc from the graph.
reverse.arc = function(x, from, to, check.cycles = TRUE, debug = FALSE) {

  arc.operations(x = x, from = from, to = to, op = "reverse",
    check.cycles = check.cycles, debug = debug)

}#REVERSE.ARC

arc.operations = function(x, from, to, op = NULL, check.cycles, update = TRUE, debug = FALSE) {

available.ops = c("set", "drop", "reverse")

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # check the op code.
  if (!(op %in% available.ops))
    stop("valid op codes are 'set', 'drop' and 'reverse'.")
  # both 'from' and 'to' are needed.
  if (missing(from) || missing(to))
    stop("'from' and/or 'to' are missing.")
  # 'from' must be a character string.
  if (!is(from, "character"))
    stop("'from' must be a character string, the label of a node.")
  # only one 'from' node is needed.
  if (length(from) > 1)
    stop("only a single 'from' node may be specified.")
  # 'to' must be a character string.
  if (!is(to, "character"))
    stop("'to' must be a character string, the label of a node.")
  # only one 'to' node is needed.
  if (length(to) > 1)
    stop("only a single 'to' node may be specified.")
  # 'from' must be a valid node label.
  if (!(from %in% names(x$nodes)))
    stop("'from' not present in the graph.")
  # 'to' must be a valid node label.
  if (!(to %in% names(x$nodes)))
    stop("'to' not present in the graph.")
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check check.cycles.
  if (!is.logical(check.cycles) || is.na(check.cycles))
    stop("check.cycles must be a logical value (TRUE/FALSE).")
  # check check.cycles.
  if (!is.logical(update) || is.na(update))
    stop("update must be a logical value (TRUE/FALSE).")

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
score = function(x, data, type = "bde", debug = FALSE) {

  available.types = c("lik", "loglik", "aic", "bic", "dir", "bde")

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # the original data set is needed.
  if (missing(data))
    stop("the data are missing.")
  # check the data for NULL/NaN/NA
  if (missing.data(data))
    stop("the data set contains NULL/NaN/NA values.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  # no score if the graph is partially directed.
  if (any(is.undirected(x$arcs)))
    stop("the graph is only partially directed.")
  # check the score label against available.types.
  if (!(type %in% available.types))
    stop("unknown score type.")
  # only discrete bayesian networks are supported.
  if (!is.data.discrete(data))
    stop("scores for continuous networks not implemented.")

  if (type %in% c("dir", "bde")) 
    return(bde.score(x = x, data = data, debug = debug))
  else if (type == "loglik")
   return(loglik.score(x = x, data = data, debug = debug))
  else if (type == "lik")
    return(exp(loglik.score(x = x, data = data, debug = debug)))
  else if (type == "aic")
    return(loglik.score(x = x, data = data, penalized = 1, debug = debug))
  else if (type == "bic")
    return(loglik.score(x = x, data = data, penalized = log(nrow(data))/2, 
      debug = debug))

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
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")

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

choose.direction = function(x, arc, data, debug = FALSE) {

  # check x's class.
  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  # check the data are there.
  if (missing(data))
    stop("the data are missing.")
  # check the data for NULL/NaN/NA.
  if (missing.data(data))
    stop("the data set contains NULL/NaN/NA values.")
  # check the arc is there.
  if (missing(arc))
    stop("the arc is missing.")
 # nodes must be a vector of character strings.
  if (!is(arc, "character"))
    stop("arc must be a vector of character strings, the labels of the nodes.")
  # exactly two nodes are needed.
  if (length(arc) != 2)
    stop("an arc if formed by exactly 2 nodes.")
  # no duplicates allowed.
  if(any(duplicated(arc)))
     stop("node labels in the arc must be unique.")
  # both elements of the arc must be valid node labels.
  if (!all(arc %in% names(x$nodes)))
    stop("node not present in the graph.")
  # check debug.
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")

  if (debug)
    cat("* testing", arc[1], "-", arc[2], "for direction.\n" )

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wonder why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.

  parents1 = parents.backend(x$arcs, arc[2], TRUE)
  a1 = conditional.test(arc[1], arc[2],
        parents1[parents1 != arc[1]],
        data = data, test = x$learning$test)

  parents2 = parents.backend(x$arcs, arc[1], TRUE)
  a2 = conditional.test(arc[2], arc[1],
        parents2[parents2 != arc[2]],
        data = data, test = x$learning$test)

  if (debug) {

    cat("  > testing", arc[1], "->", arc[2], "with conditioning set '",
      parents1[parents1 != arc[1]], "'.\n")
    cat("    > p-value is", a1, ".\n")

    cat("  > testing", arc[2], "->", arc[1], "with conditioning set '",
      parents2[parents2 != arc[2]], "'.\n")
    cat("    > p-value is", a2, ".\n")

  }#THEN

  if ((a2 < a1) && (a2 < x$learning$alpha)) {

    if (debug) cat("  @ removing", arc[1], "->", arc[2], ".\n")

    x$arcs = set.arc.direction(arc[2], arc[1], x$arcs)
    x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  }#THEN
  else if ((a1 < a2) && (a1 < x$learning$alpha)) {

    if (debug) cat("  @ removing", arc[2], "->", arc[1], ".\n")

    x$arcs = set.arc.direction(arc[1], arc[2], arcs)
    x$nodes = cache.structure(names(x$nodes), x$arcs, debug = debug)

  }#THEN
  else {

    if (debug) cat("  @ nothing to be done.\n")

  }#ELSE

  invisible(x)

}#CHOOSE.DIRECTION

