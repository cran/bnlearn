
# Parameter sanitization for arc operations.
arc.operations = function(x, from, to, op = NULL, check.cycles, update = TRUE,
    debug = FALSE) {

  available.ops = c("set", "drop", "reverse")

  # check x's class.
  check.bn(x)
  # check the op code.
  if (!(op %in% available.ops))
    stop("valid op codes are 'set', 'drop' and 'reverse'.")
  # a valid node is needed.
  check.nodes(nodes = from, graph = x, max.nodes = 1)
  # another valid node is needed.
  check.nodes(nodes = to, graph = x, max.nodes = 1)
  # 'from' must be different from 'to'.
  if (identical(from, to))
    stop("'from' and 'to' must be different from each other.")
  # check logical flags (debug, check.cycles, update).
  check.logical(debug)
  check.logical(check.cycles)
  check.logical(update)

  # add/reverse/orient the arc.
  if (op == "set") {

    if (debug)
      cat("* setting arc", from, "->", to, ".\n")

    x$arcs = set.arc.direction(from, to, x$arcs, debug = debug)

  }#THEN
  else if (op == "drop") {

    if (debug)
      cat("* dropping any arc between ", from, "and", to, ".\n")

    x$arcs = drop.arc.backend(x$arcs, c(from, to), debug = debug)

  }#THEN
  else if (op == "reverse") {

    if (debug)
      cat("* reversing any arc between ", from, "and", to, ".\n")

    x$arcs = reverse.arc.backend(from, to, x$arcs, debug = debug)

  }#THEN

  # check whether the graph is still acyclic; not needed if an arc is dropped.
  if (check.cycles && (op != "drop"))
    if (!is.acyclic(x$arcs, names(x$nodes)))
      stop("the resulting graph contains cycles.")

  # update the network structure.
  if (update) {

    # build the adjacency matrix only once.
    amat = arcs2amat(x$arcs, names(x$nodes))
    # check which nodes have to be updated.
    updated.nodes = unique(c(from, to, x$nodes[[from]]$mb, x$nodes[[to]]$mb))
    # update the chosen nodes.
    for (node in updated.nodes)
      x$nodes[[node]] = cache.partial.structure(names(x$nodes),
        target = node, amat = amat, debug = debug)

  }#THEN

  invisible(x)

}#ARC.OPERATIONS

