
# sanitize the extra arguments passed to the graph enumerations.
check.enumeration.args = function(type, N, extra) {

  # check the number of root nodes.
  if (has.argument(type, "k", enumerations.extra.args))
    check.graph.root.nodes(extra$k, N)

  # check the number of arcs.
  if (has.argument(type, "r", enumerations.extra.args))
    check.graph.narcs(extra$r, N)

  check.unused.args(extra, enumerations.extra.args[[type]])

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
