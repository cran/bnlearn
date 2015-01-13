
# return the nodes in the graph.
.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (is(x, "bn"))
    names(x$nodes)
  else
    names(x)

}#.NODES

# relabel the nodes of a graph.
.relabel = function(x, value) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid replacement node set is needed.
  nodes = .nodes(x)
  nnodes = length(nodes)
  check.nodes(nodes = value, min.nodes = nnodes, max.nodes = nnodes)

  if (is(x, "bn"))
    .relabel.bn(x, nodes, value)
  else
    .relabel.bn.fit(x, nodes, value)

}#.RELABEL

.relabel.bn = function(x, nodes, value) {

  # first, update the arcs.
  x$arcs[] = value[match(x$arcs, nodes)]
  # second, regenerate cached information.
  x$nodes = cache.structure(value, arcs = x$arcs)
  # third, the learning algorithm meta-data.
  # whitelist...
  if (!is.null(x$learning$whitelist))
    x$learning$whitelist[] = value[match(x$learning$whitelist, nodes)]
  # ... blacklist...
  if (!is.null(x$learning$blacklist))
    x$learning$blacklist[] = value[match(x$learning$blacklist, nodes)]
  # ... Castelo & Siebes prior specification...
  if (("prior" %in% names(x$learning$args)) &&
      (x$learning$args$prior == "cs")) {

    x$learning$args$beta[, c("from", "to")] =
      value[match(x$learning$args$beta[, c("from", "to")], nodes)]
    attr(x$learning$args$beta, "nodes") = value

  }#THEN
  # ... training node for classifiers...
  if (is(x, c("bn.naive", "bn.tan")))
    x$learning$args$training = value[match(x$learning$args$training, nodes)]
  # ... and the root node of TAN.
  if (is(x, "bn.tan"))
    x$learning$args$root = value[match(x$learning$args$root, nodes)]

  return(x)

}#.RELABEL.BN

.relabel.bn.fit = function(x, nodes, value) {

  for (ldist in nodes)  {

    # extract the local distribution for easy reference.
    l = x[[ldist]]
    # update the structure meta-data.
    l$node = value[match(ldist, nodes)]
    l$parents = value[match(l$parents, nodes)]
    l$children = value[match(l$children, nodes)]

    if (is(l, c("bn.fit.dnode", "bn.fit.onode"))) {

      # update the dimnames of the conditional probability table.
      if (length(l$parents) > 0)
        names(dimnames(l$prob)) = value[match(names(dimnames(l$prob)), nodes)]

    }#THEN
    else if (is(l, c("bn.fit.gnode"))){

      # update the labels of the regression coefficients.
      if (length(l$parents) > 0) {

        # take care not to mess with the intercept.
        coefs = names(l$coefficients)[names(l$coefficients) != "(Intercept)"]
        coefs = value[match(coefs, nodes)]
        names(l$coefficients)[names(l$coefficients) != "(Intercept)"] = coefs

      }#THEN

    }#ELSE

    x[ldist] = list(l)

  }#FOR

  # relabel the nodes.
  names(x) = value

  # update the training node of classifiers.
  if (is(x, c("bn.naive", "bn.tan")))
    attr(x, "training") = value[match(attr(x, "training"), nodes)]

  return(x)

}#.RELABEL.BN.FIT

# get the degree of a node.
.degree = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
    length(x$nodes[[node]]$nbr)
  else
    length(x[[node]]$parents) + length(x[[node]]$children)

}#.DEGREE


