
# this is to keep the old S3 behaviour inside the NAMESPACE.
is = function(x, class) {

    any(class(x) %in% class)

}#IS

# make bnlearn's classes known to S4.
setClass("bn")
setClass("bn.fit")
setClass("bn.naive")
setClass("bn.tan")

# return the nodes in the graph.
.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (inherits(x, "bn"))
    names(x$nodes)
  else
    names(x)

}#.NODES

setMethod("nodes", "bn", function(object) .nodes(object))
setMethod("nodes", "bn.fit", function(object) .nodes(object))
setMethod("nodes", "bn.naive", function(object) .nodes(object))
setMethod("nodes", "bn.tan", function(object) .nodes(object))

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

setMethod("degree", "bn", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.fit", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.naive", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.tan", function(object, Nodes) .degree(object, Nodes))

# convert a bn or bn.fit object into a graphNEL one.
setAs(from = "bn", to = "graphNEL", function(from, to) as.graphNEL.bn(from))
setAs(from = "bn.fit", to = "graphNEL", function(from, to) as.graphNEL.bn.fit(from))
setAs(from = "bn.naive", to = "graphNEL", function(from, to) as.graphNEL(from))
setAs(from = "bn.tan", to = "graphNEL", function(from, to) as.graphNEL(from))

# convert a bn or bn.fit object into a graphAM one.
setAs(from = "bn", to = "graphAM", function(from, to) as.graphAM.bn(from))
setAs(from = "bn.fit", to = "graphAM", function(from, to) as.graphAM.bn.fit(from))
setAs(from = "bn.naive", to = "graphAM", function(from, to) as.graphAM(from))
setAs(from = "bn.tan", to = "graphAM", function(from, to) as.graphAM(from))

# covert graph objects back to bn objects.
setAs(from = "graphNEL", to = "bn", function(from, to) as.bn.graphNEL(from))
setAs(from = "graphAM", to = "bn", function(from, to) as.bn.graphAM(from))

