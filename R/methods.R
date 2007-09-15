# compare two graphs
compare = function (r1, r2, debug = FALSE) {

  result = TRUE

  # check both objects' class.
  if (!is(r1, "bn") || !is(r2, "bn"))
    stop("both r1 and r2 must be objects of class 'bn'.")

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

        return(FALSE)

      }#THEN

      # ... and the neighbourhood.
      if (!identical(sort(r1.node$nbr), sort(r2.node$nbr))) {

        if (debug) {

          cat("* nodes in the neighbourhood of", node, "in r1 not present in r2:\n")
          print(r1.node$nbr[!(r1.node$nbr %in% r2.node$nbr)])
          cat("* nodes in the neighbourhood of", node, "in r2 not present in r1:\n")
          print(r2.node$nbr[!(r2.node$nbr %in% r1.node$nbr)])

        }#THEN

        return(FALSE)

      }#THEN

      return(TRUE)

    }

  )

  if (!all(check)) result = FALSE

  # check the arcs.
  # build two arrays of labels for easy processing.
  r1.arcs = apply(r1$arcs, 1, paste, collapse = " -> ")
  r2.arcs = apply(r2$arcs, 1, paste, collapse = " -> ")

  if (!identical(sort(r1.arcs), sort(r2.arcs))) {

    if (debug) {

      cat("* arcs in r1 not present in r2:\n")
      print(r1.arcs[!(r1.arcs %in% r2.arcs)])
      cat("* arcs in r2 not present in r1:\n")
      print(r2.arcs[!(r2.arcs %in% r1.arcs)])

    }#THEN

    result = FALSE

  }#THEN

  result

}#COMPARE

# return the markov blanket of a node
mb = function(x, node) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")

  x$nodes[[node]]$mb

}#MB

# return the neighbourhood of a node
nbr = function(x, node) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")

  x$nodes[[node]]$nbr

}#NBR

# return the arcs in the graph
arcs = function(x) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  x$arcs

}#ARCS

# return the nodes in the graph
nodes = function(x) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  names(x$nodes)

}#NODES

# build an adjacency matrix from a graph.
amat = function(x) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")

  arcs2amat(x$arcs, names(x$nodes))

}#AMAT

# get the parents of a node
parents = function(x, node) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  if (!(node %in% names(x$nodes)))
    stop("node not present in the graph.")

   parents.backend(x$arcs, node)

}#PARENTS

choose.direction = function(x, arc, data, debug = FALSE) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  if (!all(arc %in% names(x$nodes)))
    stop("node not present in the graph.")
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")

  if (debug)
    cat("* testing", arc[1], "-", arc[2], "for direction.\n" )

  # you can't help but notice nodes connected by undirected arcs are
  # included, too? wondwer why?
  # because if they, too, are parents of the node to be tested
  # they _do_ belong there; if they are not, the node distribution
  # does not depend on them so they are largely irrelevant.

  parents1 = parents.backend(x$arcs, arc[2], TRUE)
  a1 = conditional.test(arc[1], arc[2], 
        parents1[parents1 != arc[1]], 
        data = data, test = x$test)

  parents2 = parents.backend(x$arcs, arc[1], TRUE)
  a2 = conditional.test(arc[2], arc[1], 
        parents2[parents2 != arc[2]], 
        data = data, test = x$test)

  if (debug) {

    cat("  > testing", arc[1], "->", arc[2], "with conditioning set '", 
      parents1[parents1 != arc[1]], "'.\n")
    cat("    > p-value is", a1, ".\n")

    cat("  > testing", arc[2], "->", arc[1], "with conditioning set '", 
      parents2[parents2 != arc[2]], "'.\n")
    cat("    > p-value is", a2, ".\n")

  }#THEN

  if (a2 < a1) {

    if(debug) cat("  @ removing", arc[1], "->", arc[2], ".\n")

    x$arcs = x$arcs[!is.row.equal(x$arcs,arc), ]

  }#THEN
  else if (a1 < a2) {

    if(debug) cat("  @ removing", arc[2], "->", arc[1], ".\n")

    x$arcs = x$arcs[!is.row.equal(x$arcs, arc[c(2,1)]), ]

  }#ELSE

  invisible(x)

}#CHOOSE.DIRECTION

loglik = function(x, data, debug = FALSE) {

  if (!is(x, "bn"))
    stop("x must be an object of class 'bn'.")
  if (!is.logical(debug) || is.na(debug))
    stop("debug must be a logical value (TRUE/FALSE).")
  if (!x$discrete)
    stop("continuous networks are not supported at present.")

  if (any(which.undirected(x$arcs)))
    stop("the graph is only partially directed.")

  sum(sapply(names(x$nodes), loglik.node, x = x, data = data, debug = debug))

}#LOGLIK

