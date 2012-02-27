
# make bnlearn's classes known to S4.
setClass("bn")
setClass("bn.fit")

# return the nodes in the graph.
setMethod("nodes", "bn", function(object) names(object$nodes))
setMethod("nodes", "bn.fit", function(object) names(object))

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

}#DEGREE

setMethod("degree", "bn", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.fit", function(object, Nodes) .degree(object, Nodes))

# convert a bn or bn.fit object into a graphNEL one.
setAs(from = "bn", to = "graphNEL", function(from, to) as.graphNEL.bn(from))
setAs(from = "bn.fit", to = "graphNEL", function(from, to) as.graphNEL.bn.fit(from))

# convert a bn or bn.fit object into a graphAM one.
setAs(from = "bn", to = "graphAM", function(from, to) as.graphAM.bn(from))
setAs(from = "bn.fit", to = "graphAM", function(from, to) as.graphAM.bn.fit(from))

# covert graph objects back to bn objects.
setAs(from = "graphNEL", to = "bn", function(from, to) as.bn.graphNEL(from))
setAs(from = "graphAM", to = "bn", function(from, to) as.bn.graphAM(from))

