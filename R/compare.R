
# compare two bayesian network structures.
compare.backend = function(target, current, debug = FALSE) {

  result = TRUE

  # check the two graphs have the same nodes.
  target.nodes = names(target$nodes)
  current.nodes = names(current$nodes)

  if (!identical(sort(target.nodes), sort(current.nodes))) {

    if (debug) {

      cat("* nodes in target not present in current:\n")
      print(target.nodes[!(target.nodes %in% current.nodes)])
      cat("* nodes in current not present in target:\n")
      print(current.nodes[!(current.nodes %in% target.nodes)])

    }#THEN

    return(FALSE)

  }#THEN

  # for each node check ...
  check = sapply(names(target$nodes),

    function(node) {

      node.result = TRUE
      target.node = target$nodes[[node]]
      current.node = current$nodes[[node]]

      # ... the markov blanket ...
      if (!identical(sort(target.node$mb), sort(current.node$mb))) {

        if (debug) {

          cat("* nodes in the markov blanket of", node, "in target not present in current:\n")
          print(target.node$mb[!(target.node$mb %in% current.node$mb)])
          cat("* nodes in the markov blanket of", node, "in current not present in target:\n")
          print(current.node$mb[!(current.node$mb %in% target.node$mb)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... and the neighbourhood ...
      if (!identical(sort(target.node$nbr), sort(current.node$nbr))) {

        if (debug) {

          cat("* nodes in the neighbourhood of", node, "in target not present in current:\n")
          print(target.node$nbr[!(target.node$nbr %in% current.node$nbr)])
          cat("* nodes in the neighbourhood of", node, "in current not present in target:\n")
          print(current.node$nbr[!(current.node$nbr %in% target.node$nbr)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... the parents ...
      if (!identical(sort(target.node$parents), sort(current.node$parents))) {

        if (debug) {

          cat("* parents of", node, "in target not present in current:\n")
          print(target.node$parents[!(target.node$parents %in% current.node$parents)])
          cat("* parents of", node, "in current not present in target:\n")
          print(current.node$parents[!(current.node$parents %in% target.node$parents)])

        }#THEN

        node.result = FALSE

      }#THEN

      # ... and the children.
      if (!identical(sort(target.node$children), sort(current.node$children))) {

        if (debug) {

          cat("* children of", node, "in target not present in current:\n")
          print(target.node$children[!(target.node$children %in% current.node$children)])
          cat("* children of", node, "in current not present in target:\n")
          print(current.node$children[!(current.node$children %in% target.node$children)])

        }#THEN

        node.result = FALSE

      }#THEN

      return(node.result)

    }

  )

  if (!all(check)) result = FALSE

  # check directed arcs.
  target.arcs = apply(target$arcs[which.directed(target$arcs), , drop = FALSE], 1, paste, collapse = " -> ")
  current.arcs = apply(current$arcs[which.directed(current$arcs), , drop = FALSE], 1, paste, collapse = " -> ")

  if (!identical(sort(target.arcs), sort(current.arcs))) {

    if (debug) {

      cat("* directed arcs in target not present in current:\n")
      print(target.arcs[!(target.arcs %in% current.arcs)])
      cat("* directed arcs in current not present in target:\n")
      print(current.arcs[!(current.arcs %in% target.arcs)])

    }#THEN

    result = FALSE

  }#THEN

  # check undirected arcs.
  target.arcs = apply(target$arcs[which.undirected(target$arcs), , drop = FALSE], 1, paste, collapse = " - ")
  current.arcs = apply(current$arcs[which.undirected(current$arcs), , drop = FALSE], 1, paste, collapse = " - ")

  if (!identical(sort(target.arcs), sort(current.arcs))) {

    if (debug) {

      cat("* undirected arcs in target not present in current:\n")
      print(target.arcs[!(target.arcs %in% current.arcs)])
      cat("* undirected arcs in current not present in target:\n")
      print(current.arcs[!(current.arcs %in% target.arcs)])

    }#THEN

    result = FALSE

  }#THEN

  result

}#COMPARE

