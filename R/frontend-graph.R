
# get the root nodes of a graph.
root.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  root.leaf.nodes(x, leaf = FALSE)

}#ROOT.NODES

# get the leaf nodes of a graph.
leaf.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  root.leaf.nodes(x, leaf = TRUE)

}#LEAF.NODES

# check if a graph is acyclic.
acyclic = function(x, directed, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)

  # fitted bayesian networks are always acylic.
  if (class(x) == "bn.fit")
    return(TRUE)

  if (missing(directed)) {

    is.acyclic(x$arcs, names(x$nodes), debug = debug)

  }#THEN
  else {

    # check directed.
    check.logical(directed)

    is.acyclic.backend(arcs = x$arcs, nodes = names(x$nodes),
      directed = directed, debug = debug)

  }#ELSE

}#ACYCLIC

# check if a graph is directed.
directed = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  # fitted bayesian networks are always directed.
  if (class(x) == "bn.fit")
    return(TRUE)

  is.dag(x$arcs, names(x$nodes))

}#DIRECTED

# check if there's a path between two specific nodes.
path = function(x, from, to, direct = TRUE, underlying.graph = FALSE,
    debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = from, graph = x, max.nodes = 1)
  # another valid node is needed.
  check.nodes(nodes = to, graph = x, max.nodes = 1)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check underlying.path.
  check.logical(underlying.graph)
  # check debug.
  check.logical(debug)

  if (class(x) == "bn") {

    nodes = names(x$nodes)
    amat = arcs2amat(x$arcs, nodes)

  }#THEN
  else {

    nodes = names(x)
    amat = arcs2amat(fit2arcs(x), nodes)

  }#ELSE

  has.path(from, to, nodes = nodes, amat = amat, exclude.direct = !direct,
    underlying.graph = underlying.graph, debug = debug)

}#PATH

# return the partial node ordering implied by the graph structure.
node.ordering = function(x, debug = FALSE) {

  # check x's class.
  check.bn.or.fit(x)
  # check debug.
  check.logical(debug)
  # no model string if the graph is partially directed.
  if (class(x) == "bn")
    if (is.pdag(x$arcs, names(x$nodes)))
      stop("the graph is only partially directed.")

  schedule(x, debug = debug)

}#NODE.ORDERING

# generate a valid blacklist from a partial node ordering.
ordering2blacklist = function(nodes) {

  # check the node labels.
  check.nodes(nodes, min.nodes = 3)

  do.call(rbind,

    sapply(seq(from = 1, to = length(nodes)),
      function(i) {

        cbind(from = rep(nodes[i], i - 1), to = nodes[(1:i) - 1])

    })
  )

}#ORDERING2BLACKLIST

# return the skeleton of a (partially) directed graph
skeleton = function(x) {

  # check x's class.
  check.bn(x)

  dag2ug.backend(x, moral = FALSE)

}#SKELETON

# return a complete orientation of a graph.
pdag2dag = function(x, ordering) {

  # check x's class.
  check.bn(x)
  # check the node ordering.
  check.nodes(ordering, graph = x, min.nodes = length(x$nodes),
    max.nodes = length(x$nodes))

  pdag2dag.backend(x, ordering)

}#PDAG2DAG

# compare two bayesian network structures.
compare = function (r1, r2, debug = FALSE) {

  result = TRUE

  # check both objects' class.
  check.bn(r1)
  check.bn(r2)
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
  r1.arcs = apply(r1$arcs[which.directed(r1$arcs), , drop = FALSE], 1, paste, collapse = " -> ")
  r2.arcs = apply(r2$arcs[which.directed(r2$arcs), , drop = FALSE], 1, paste, collapse = " -> ")

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
  r1.arcs = apply(r1$arcs[which.undirected(r1$arcs), , drop = FALSE], 1, paste, collapse = " - ")
  r2.arcs = apply(r2$arcs[which.undirected(r2$arcs), , drop = FALSE], 1, paste, collapse = " - ")

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

