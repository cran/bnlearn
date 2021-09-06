
# return the markov blanket of a node.
mb = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
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

  if (is(x, "bn"))
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

  if (is(x, "bn"))
    x$nodes[[node]]$parents
  else
    x[[node]]$parents

}#PARENTS

# add one or more parents to a node.
"parents<-" = function(x, node, debug = FALSE, value) {

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
  to.be.dropped = x$nodes[[node]]$parents[x$nodes[[node]]$parents %!in% value]
  # add only the nodes that were not already there.
  to.be.added = value[value %!in% x$nodes[[node]]$parents]

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
          check.cycles = FALSE, check.illegal = FALSE, update = FALSE,
          debug = debug)

  }#FOR

  # adding!
  for (p in to.be.added) {

    x = arc.operations(x = x, from = p, to = node, op = "set",
          check.cycles = TRUE, check.illegal = TRUE, update = FALSE,
          debug = debug)

  }#FOR

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs, debug = debug)

  return(x)

}#PARENTS<-

# get the children of a node.
children = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
    x$nodes[[node]]$children
  else
    x[[node]]$children

}#CHILDREN

# add one or more children to a node.
"children<-" = function(x, node, debug = FALSE, value) {

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
  to.be.dropped = x$nodes[[node]]$children[x$nodes[[node]]$children %!in% value]
  # add only the nodes that were not already there.
  to.be.added = value[value %!in% x$nodes[[node]]$children]

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
          check.cycles = FALSE, check.illegal = FALSE, update = FALSE,
          debug = debug)

  }#FOR

  # adding!
  for (child in to.be.added) {

    x = arc.operations(x = x, from = node, to = child, op = "set",
          check.cycles = TRUE, check.illegal = TRUE, update = FALSE,
          debug = debug)

  }#FOR

  # update the network structure.
  x$nodes = cache.structure(names(x$nodes), arcs = x$arcs, debug = debug)

  return(x)

}#CHILDREN<-

# get the spouses of a node.
spouses = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  chld = children(x, node)
  sp = unique(unlist(lapply(chld, parents, x = x)))

  return(setdiff(as.character(sp), node))

}#SPOUSES

# get the ancestors of a node.
ancestors = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  # the first element is the node itself, which is not its own ancestor.
  return(topological.ordering(x, start = node, reverse = TRUE)[-1])

}#ANCESTORS

# get the descendants of a node.
descendants = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  # the first element is the node itself, which is not its own ancestor.
  return(topological.ordering(x, start = node)[-1])

}#DESCENDANTS

# get the in-degree of a node.
in.degree = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
    length(x$nodes[[node]]$parents)
  else
    length(x[[node]]$parents)

}#IN.DEGREE

# get the out-degree of a node.
out.degree = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
    length(x$nodes[[node]]$children)
  else
    length(x[[node]]$children)

}#OUT.DEGREE

# add a node to a network structure.
add.node = function(x, node) {

  # check x's class.
  check.bn(x)
  # a valid node (that is not already in the graph) is needed.
  if (!is.string(node))
    stop("'node' should be a single character string.")
  if (node %in% names(x$nodes))
    stop("node ", node, " is already present in the graph.")

  # add the node, without adding arcs.
  x$nodes[[node]] = list(mb = character(0), nbr = character(0),
                         parents = character(0), children = character(0))

  return(x)

}#ADD.NODE

# remove a node from a network structure.
remove.node = function(x, node) {

  # check x's class.
  check.bn(x)
  # the resulting graphs should have at least one node.
  if (length(x$nodes) == 1)
    stop("trying to remove the only node in the graph.")
  # a valid node (that is already in the graph) is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  # removing a single node is a special case of creating a subgraph.
  subgraph.backend(x = x, nodes = setdiff(names(x$nodes), node),
    preserve.learning = TRUE)

}#REMOVE.NODE

# rename the nodes of a network structure.
rename.nodes = function(x, names) {

  .relabel(x, value = names)

}#RENAME.NODES
