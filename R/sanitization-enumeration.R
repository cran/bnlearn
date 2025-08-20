
# sanitize the extra arguments passed to the graph enumerations.
check.enumeration.args = function(type, extra) {

  # check the number of nodes.
  if (has.argument(type, "nodes", enumerations.extra.args)) {

    if (!is.positive.vector(extra[["nodes"]]))
      stop("'nodes' must be positive integers, the number(s) of nodes in the graph.")

  }#THEN

  # check the number of root nodes.
  if (has.argument(type, "k", enumerations.extra.args))
    check.graph.root.nodes(extra[["k"]], N = extra[["nodes"]])

  # check the number of arcs.
  if (has.argument(type, "r", enumerations.extra.args))
    check.graph.narcs(extra[["r"]], N = extra[["nodes"]])

  # check the markov equivalence class.
  if (has.argument(type, "eqclass", enumerations.extra.args))
    valid.cpdag.backend(extra[["eqclass"]])

  # warn about and remove unused arguments.
  extra = check.unused.args(extra, enumerations.extra.args[[type]])

  return(extra)

}#CHECK.ENUMERATION.ARGS

# check the number of root nodes.
check.graph.root.nodes = function(k, N) {

  if (is.null(k))
    stop("unspecified number of root nodes.")

  if (!is.positive.integer(k))
    stop("the number of root nodes must be a positive integer.")
  if (k > max(N))
    warning("the number of root nodes is larger than the number of nodes.")

}#CHECK.GRAPH.ROOT.NODES

# check the number of arcs.
check.graph.narcs = function(r, N) {

  if (is.null(r))
    stop("unspecified number of arcs.")

  if (!is.non.negative.integer(r))
    stop("the number of arcs must be a positive integer.")
  if (r > choose(max(N), 2))
    warning("the number of arcs is larger than the maximum possible number of arcs.")

}#CHECK.GRAPH.NARCS
