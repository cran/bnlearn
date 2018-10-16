
# check nodes (not necessarily from a bn object).
check.nodes = function(nodes, graph = NULL, min.nodes = 1, max.nodes = Inf) {

  # a node is needed.
  if (missing(nodes))
    stop("no node specified.")
  # nodes must be a vector of character strings.
  if (!is.string.vector(nodes))
    stop("nodes must be a vector of character strings, the labels of the nodes.")
  # no duplicates allowed.
  if (any(duplicated(nodes)))
     stop("node labels must be unique.")
  # required maximum number of nodes.
  if (length(nodes) > max.nodes)
    stop("at most ", max.nodes, " node(s) needed.")
  # required minimum number of nodes (usually 1).
  if (length(nodes) < min.nodes)
    stop("at least ", min.nodes, " node(s) needed.")
  # node must be a valid node label.
  if (!is.null(graph)) {

    if (is(graph, "bn"))
      invalid.nodes = nodes %!in% names(graph$nodes)
    else if (is(graph, "bn.fit"))
      invalid.nodes = nodes %!in% names(graph)
    else if (is.character(graph))
      invalid.nodes = nodes %!in% graph

    if (any(invalid.nodes))
        stop("invalid node(s)",
          paste0(" '", nodes[invalid.nodes], "'"), ".")

  }#THEN

}#CHECK.NODES

# check an arc set.
check.arcs = function(arcs, nodes) {

  # sanitize the set of arcs.
  if (is(arcs, "matrix") || is(arcs, "data.frame")) {

     if (dim(arcs)[2] != 2)
       stop("arc sets must have two columns.")
     if (!all(sapply(arcs, class) == "character"))
       stop("node labels in arc sets must be character strings.")

     if (is.data.frame(arcs))
       arcs = as.matrix(cbind(as.character(arcs[, 1]),
         as.character(arcs[, 2])))

     # be sure to set the column names.
     dimnames(arcs) = list(c(), c("from", "to"))

  }#THEN
  else if (is.character(arcs)) {

    # if there is an even number of labels fit them into a 2-column matrix.
    if ((length(arcs) %% 2) != 0)
      stop("arc sets must have two columns.")

    arcs = matrix(arcs, ncol = 2, byrow = TRUE,
              dimnames = list(c(), c("from", "to")))

  }#THEN
  else {

     stop("an arc set must be a matrix or data.frame with two columns.")

  }#ELSE

  # nodes must be valid node labels.
  if (any(arcs %!in% nodes))
    stop("node(s)", paste0(" '", unique(arcs[arcs %!in% nodes]), "'"),
         " not present in the graph.")

  # remove duplicate arcs.
  arcs = unique.arcs(arcs, nodes, warn = TRUE)

  # check there are no loops among the arcs.
  loop = (arcs[, "from"] == arcs[, "to"])

  if (any(loop))
    stop("invalid arcs that are actually loops:\n",
         paste("  ", arcs[loop, 1], "->", arcs[loop, 2], "\n"))

  return(arcs)

}#CHECK.ARCS

# check node groups and partitions.
check.node.groups = function(groups, graph = NULL) {

  # check the node labels.
  if (is.list(groups)) {

    vname = deparse(substitute(groups))

    if (length(groups) == 0)
      stop("no ", vname, " specified.")
    if (!all(sapply(groups, is.character)))
      stop("node labels in '", vname, "' must be character strings.")
    if (!all(sapply(groups, length) > 0))
      stop("all ", vname, " must contain at least one node label.")

    check.nodes(unlist(groups), graph = graph)

  }#THEN
  else {

    check.nodes(groups, graph = graph)

  }#ELSE

}#CHECK.NODE.GROUPS
