
## generics for the S3 methods below (and others in the frontend-* files).

sigma = function(object, ...) {

  UseMethod("sigma")

}#SIGMA

as.grain = function(x) {

  UseMethod("as.grain", x)

}#AS.GRAIN

as.bn = function(x) {

  UseMethod("as.bn", x)

}#AS.BN

as.bn.fit = function(x) {

  UseMethod("as.bn.fit", x)

}#AS.BN.FIT

as.graphNEL = function(x) {

  UseMethod("as.graphNEL", x)

}#AS.GRAPHNEL

as.graphAM = function(x) {

  UseMethod("as.graphAM", x)

}#AS.GRAPHAM

# convert a "bn.fit" object from bnlearn into a "grain" object.
as.grain.bn.fit = function(x) {

  # check whether gRain is loaded.
  if (!requireNamespace("gRain"))
    stop("this function requires the gRain package.")

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

      values[is.na(values)] <- 1/dim(values)[1]

    }#THEN

    cpt[[node]] = gRain::cptable(formula(f), values = values, levels = levels)

  }#FOR

  return(gRain::grain.CPTspec(gRain::compileCPT(cpt)))

}#AS.GRAIN.BN.FIT

# convert a "grain" object from gRain into a "bn.fit" object.
as.bn.fit.grain = function(x) {

  # check whether gRain is loaded.
  if (!requireNamespace("gRain"))
    stop("this function requires the gRain package.")

  if (!is(x, "grain"))
    stop("x must be an object of class 'grain'.")

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

  return(structure(fitted, class = c("bn.fit", guess.fitted.class(fitted))))

}#AS.BN.FIT.GRAIN

# convert a "grain" object from gRain into a "bn" object.
as.bn.grain = function(x) {

  # check whether gRain is loaded.
  if (!requireNamespace("gRain"))
    stop("this function requires the gRain package.")

  if (!is(x, "grain"))
    stop("x must be an object of class 'grain'.")

  return(bn.net(as.bn.fit(x)))

}#AS.BN.GRAIN

# convert a "bn" or "bn.fit" object into a "graphNEL" object.
as.graphNEL.bn = function(x) {

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
as.bn.graphAM = function(x) {

  # check whether graph is loaded.
  if (!requireNamespace("graph"))
    stop("this function requires the graph package.")

  nodes = x@nodes
  arcs = amat2arcs(x@adjMat, nodes)

  res = empty.graph(nodes)
  res$arcs = arcs
  res$nodes = cache.structure(nodes, arcs)

  return(res)

}#AS.BN.GRAPHAM

# convert a "graphNEL" object into a "bn" object.
as.bn.graphNEL = function(x) {

  # check whether graph is loaded.
  if (!requireNamespace("graph"))
    stop("this function requires the graph package.")

  nodes = x@nodes
  arcs = elist2arcs(graph::edges(x))

  res = empty.graph(nodes)
  res$arcs = arcs
  res$nodes = cache.structure(nodes, arcs)

  return(res)

}#AS.BN.GRAPHNEL

