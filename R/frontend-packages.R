
## generics for the S3 methods below (and others in the frontend-* files).

sigma = function(object, ...) {

  UseMethod("sigma")

}#SIGMA

as.grain = function(x) {

  UseMethod("as.grain")

}#AS.GRAIN

as.bn = function(x, ...) {

  UseMethod("as.bn")

}#AS.BN

as.bn.fit = function(x, ...) {

  UseMethod("as.bn.fit")

}#AS.BN.FIT

as.graphNEL = function(x) {

  UseMethod("as.graphNEL")

}#AS.GRAPHNEL

as.graphAM = function(x) {

  UseMethod("as.graphAM")

}#AS.GRAPHAM

as.prediction = function(x, ...) {

  UseMethod("as.prediction")

}#AS.PREDICTION

as.lm = function(x, ...) {

  UseMethod("as.lm")

}#AS.LM

# convert a "bn.fit" object from bnlearn into a "grain" object.
as.grain.bn.fit = function(x) {

  # check whether gRain is loaded.
  check.and.load.package("gRain")

  # check x's class.
  check.fit(x)
  # check whether x is a discrete fitted network.
  if (is(x, c("bn.fit.onet", "bn.fit.donet")))
    warning("the gRain package does not support ordinal networks, disregarding the ordering of the levels.")
  else if (!is(x, "bn.fit.dnet"))
    stop("the gRain package only supports discrete networks.")

  cpt = vector(length(x), mode = "list")
  names(cpt) = names(x)

  for (node in names(x)) {

    parents = x[[node]][["parents"]]

    if (length(parents) == 0)
      f = paste("~", node)
    else
      f = paste("~", paste(c(node, parents), collapse = "+"))

    values = x[[node]][["prob"]]
    levels = dimnames(values)[[1]]

    # gRain requires completely specified CPTs, while bnlearn does not.
    if (any(is.na(values))) {

      warning("NaN conditional probabilities in ", node,
        ", replaced with a uniform distribution.")

      values[is.na(values)] = 1/dim(values)[1]

    }#THEN

    cpt[[node]] = gRain::cptable(formula(f), values = values, levels = levels)

  }#FOR

  return(gRain::grain.CPTspec(gRain::compileCPT(cpt)))

}#AS.GRAIN.BN.FIT

# convert a "grain" object from gRain into a "bn.fit" object.
as.bn.fit.grain = function(x, ...) {

  # check whether gRain is loaded.
  check.and.load.package("gRain")

  if (!is(x, "grain"))
    stop("x must be an object of class 'grain'.")
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  grain.get.children = function(node) {

    names(which(sapply(fitted, function(x) { node %in% x$parents } )))

  }#GRAIN.GET.CHILDREN

  nodes = names(x$cptlist)

  fitted = vector(length(nodes), mode = "list")
  names(fitted) = nodes

  # first pass: get parents and CPTs.
  for (node in nodes) {

    prob = structure(as.table(x[["cptlist"]][[node]]), class = "table")
    parents = names(dimnames(prob))[-1]

    # marginal tables have no dimension names in bnlearn.
    prob = cptattr(prob)

    fitted[[node]] = structure(list(node = node, parents = parents,
                       children = NULL, prob = prob), class = "bn.fit.dnode")

  }#FOR

  # second pass: get the children.
  for (node in nodes) {

    fitted[[node]][["children"]] = grain.get.children(node)

  }#FOR

  return(structure(fitted, class = c("bn.fit", determine.fitted.class(fitted))))

}#AS.BN.FIT.GRAIN

# convert a "grain" object from gRain into a "bn" object.
as.bn.grain = function(x, ..., check.cycles = TRUE) {

  # check whether gRain is loaded.
  check.and.load.package("gRain")

  if (!is(x, "grain"))
    stop("x must be an object of class 'grain'.")
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  res = bn.net(as.bn.fit.grain(x))

  # check whether the the graph contains directed cycles.
  if (check.cycles)
    if (!is.acyclic(nodes = names(res$nodes), arcs = res$arcs, debug = FALSE,
          directed = TRUE))
      stop("the grain object contains directed cycles.")

  return(res)

}#AS.BN.GRAIN

# convert a "bn" or "bn.fit" object into a "graphNEL" object.
as.graphNEL.bn = function(x) {

  # check whether graph is loaded.
  check.and.load.package("graph")

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn")) {

    nodes = names(x$nodes)
    arcs = x$arcs

  }#THEN
  else{

    nodes = names(x)
    arcs = fit2arcs(x)

  }#ELSE

  new("graphNEL", nodes = nodes, edgeL = arcs2elist(arcs, nodes),
                  edgemode = 'directed')

}#AS.GRAPHNEL.BN

as.graphNEL.bn.fit = as.graphNEL.bn

# convert a "bn" or "bn.fit" object into a "graphAM" object.
as.graphAM.bn = function(x) {

  # check whether graph is loaded.
  if (!requireNamespace("graph"))
    stop("this function requires the graph package.")

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn")) {

    nodes = names(x$nodes)
    arcs = x$arcs

  }#THEN
  else{

    nodes = names(x)
    arcs = fit2arcs(x)

  }#ELSE

  new("graphAM", adjMat = arcs2amat(arcs, nodes), edgemode = 'directed')

}#AS.GRAPHAM.BN

as.graphAM.bn.fit = as.graphAM.bn

# convert a "graphAM" object into a "bn" object.
as.bn.graphAM = function(x, ..., check.cycles = TRUE) {

  # check whether graph is loaded.
  check.and.load.package("graph")
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  adjMat = x@adjMat
  if (is(adjMat, "double"))
    storage.mode(adjMat) = "integer"
  nodes = colnames(adjMat)
  arcs = amat2arcs(adjMat, nodes)

  # check whether the the graph contains directed cycles.
  if (check.cycles)
    if (!is.acyclic(nodes = nodes, arcs = arcs, debug = FALSE, directed = TRUE))
      stop("the graphAM object contains directed cycles.")

  res = empty.graph(nodes)
  res$arcs = arcs
  res$nodes = cache.structure(nodes, arcs)

  return(res)

}#AS.BN.GRAPHAM

# convert a "graphNEL" object into a "bn" object.
as.bn.graphNEL = function(x, ..., check.cycles = TRUE) {

  # check whether graph is loaded.
  check.and.load.package("graph")
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  nodes = x@nodes
  arcs = elist2arcs(graph::edges(x))

  # check whether the the graph contains directed cycles.
  if (check.cycles)
    if (!is.acyclic(nodes = nodes, arcs = arcs, debug = FALSE, directed = TRUE))
      stop("the graphNEL object contains directed cycles.")

  res = empty.graph(nodes)
  res$arcs = arcs
  res$nodes = cache.structure(nodes, arcs)

  return(res)

}#AS.BN.GRAPHNEL

# convert a "pcAlgo" object into a "bn" object.
as.bn.pcAlgo = function(x, ..., check.cycles = TRUE) {

  as.bn.graphNEL(x@graph, ..., check.cycles = check.cycles)

}#AS.BN.PCALGO

# generate the input for a ROC curve from arc strengths.
as.prediction.bn.strength = function(x, true, ..., consider.direction = TRUE) {

  # check whether ROCR is loaded.
  check.and.load.package("ROCR")

  # check whether true is a bn object.
  check.bn(true)
  # check whether it encodes a completely directed graph.
  nodes = names(true$nodes)
  if (is.pdag(true$arcs, nodes))
    stop("the graph is only partially directed.")
  # check the arc strengths.
  check.bn.strength(x, nodes = nodes, valid = c("bootstrap", "bayes-factor"))

  # check consider.direction.
  check.logical(consider.direction)
   # warn about unused arguments.
  check.unused.args(list(...), character(0))

  # create a data frame with all possible arcs, and merge the arc strengths
  # (with or without taking the direction into account.)
  dd = structure(data.frame(subsets(nodes, 2)), names = c("from", "to"))

  if(consider.direction) {

    true.arcs.mat = true$arcs
    x$pred = x$strength * x$direction

  }#THEN
  else {

    true.arcs.mat = arcs(skeleton(true))
    x$pred = x$strength

  }#ELSE

  # add a 0/1 indicator that encodes the arcs in the true network.
  dd = merge(dd, x[, c("from", "to", "pred")], by = c("from", "to"))
  true.arcs = data.frame(true.arcs.mat, label = 1)
  dd = merge(dd, true.arcs, by = c("from", "to"), all.x = TRUE)
  dd$label[is.na(dd$label)] = 0

  # export the (true arcs, arc probabilities) pairs to the ROCR package.
  pred = ROCR::prediction(dd$pred, dd$label)

  return(pred)

}#AS.PREDICTION.BN.STRENGTH

# fit the local distributions of a (Gaussian) network using lm().
as.lm.bn = function(x, data, ...) {

  # check x's class.
  check.bn.or.fit(x)
  # check the data.
  check.data(data, allow.missing = TRUE, stop.if.all.missing = TRUE)

  if (is(x, "bn")) {

    # check whether the data agree with the bayesian network.
    check.bn.vs.data(x, data)
    # no parameters if the network structure is only partially directed.
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

    nodes = names(x$nodes)

  }#THEN
  else {

    # check whether the data agree with the bayesian network.
    check.fit.vs.data(x, data)
    # only Gaussian networks are supported.
    if (!is(x, "bn.fit.gnet"))
      stop("only Gaussian networks are supported.")

    nodes = names(x)

  }#ELSE

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  smartSapply(NULL, nodes, function(n) {

    if (is(x, "bn"))
      node = x$nodes[[n]]
    else
      node = x[[n]]

    node$node = n
    lm.refit.node(node, data)

  })

}#AS.LM.BN

# refit the local distributions of a Gaussian network using lm().
as.lm.bn.fit = as.lm.bn

# refit a single local distribution from a Gaussian network using lm().
as.lm.bn.fit.gnode = function(x, data, ...) {

  # check the data.
  check.data(data, allow.missing = TRUE, stop.if.all.missing = TRUE)
  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  lm.refit.node(x, data)

}#AS.LM.BN.FIT.GNODE

