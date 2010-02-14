
# return the nodes in the graph.
nodes = function(x) {

  # check x's class (beware of graphviz).
  if (class(x) %in% c("graphAM", "graphNEL"))
    return(graph:::nodes(x))
  else
    check.bn.or.fit(x)

  if (class(x) == "bn")
    names(x$nodes)
  else
    names(x)

}#NODES

# return the markov blanket of a node.
mb = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (class(x) == "bn")
    x$nodes[[node]]$mb
  else
    mb.fitted(x, node)

}#MB

# return the neighbourhood of a node.
nbr = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (class(x) == "bn")
    x$nodes[[node]]$nbr
  else
    unique(c(x[[node]]$parents, x[[node]]$children))

}#NBR

# get the parents of a node.
parents = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (class(x) == "bn")
    x$nodes[[node]]$parents
  else
    x[[node]]$parents

}#PARENTS

# add one or more parents to a node.
"parents<-" <- function(x, node, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)
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
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs, debug = debug)

  x

}#PARENTS<-

# get the children of a node.
children = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (class(x) == "bn")
    x$nodes[[node]]$children
  else
    x[[node]]$children

}#CHILDREN

# add one or more children to a node.
"children<-" <- function(x, node, debug = FALSE, value) {

  # check x's class.
  check.bn(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)
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
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs, debug = debug)

  x

}#CHILDREN<-

