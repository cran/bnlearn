# check a matrix for symmetry.
is.symmetric = function(m) {

  # kill all the names; identical may return false otherwise.
  colnames(m) = rownames(m) = NULL

  identical(m, t(m))

}#IS.SYMMETRIC

# check whether a graph is completely directed (has no undirected arcs).
is.dag = function(arcs) {

  ifelse(nrow(arcs) > 0, !any(is.undirected(arcs)), TRUE)

}#IS.DAG

# alias of is.dag for readable code.
is.pdag = function(arcs) !is.dag(arcs)

# generic is.acyclic backend (calls the right one for the graph at hand).
is.acyclic = function(arcs, nodes, debug = FALSE) {

  if (any(is.undirected(arcs)))
    is.pdag.acyclic(arcs = arcs, nodes = nodes, debug = debug)
  else
    is.dag.acyclic(arcs = arcs, nodes = nodes, debug = debug)

}#IS.ACYCLIC

# check whether a (directed) graph is acyclic.
is.dag.acyclic = function(arcs, nodes, debug = FALSE) {

  # no arcs, no cycles.
  if (nrow(arcs) == 0) return(TRUE)

  amat = arcs2amat(arcs, nodes)

  if (debug)
    cat("* checking whether the (directed) graph is acyclic.\n")

  good = nodes

  # to be part of a loop, a node needs to have both a parent
  # and a children; remove (iteratively) all the nodes which
  # do not fit this description.
  for (i in seq(1, length(nodes))) {

    if (debug)
      cat("  > nodes which may be part of a cycle: '", good, "'\n")

    rs = rowSums(amat[good, good, drop = FALSE])
    cs = colSums(amat[good, good, drop = FALSE])

    # iteratively remove nodes which have either no parents
    # or no children.
    good = intersect(names(rs[rs > 0]),
                     names(cs[cs > 0]))

    # a loop requires at least three distinct nodes.
    if (length(good) < 3) {

      if (debug) cat("  @ no cycle found.\n")

      return(TRUE)

    }#THEN

  }#THEN

  if (debug)
    cat("  @ these nodes are part of one or more cycles: '", good, "'\n")

  return(FALSE)

}#IS.DIRECTED.ACYCLIC

# check whether a (partially directed) graph is acyclic.
is.pdag.acyclic = function(arcs, nodes, debug = FALSE) {

  # no arcs, no cycles.
  if (nrow(arcs) == 0) return(TRUE)

  # ignore arcs coming from root nodes or going into leaf nodes;
  # they cannot be part of a cycle beacuse of the lack of an incoming
  # or outgoing arc (respectively).
  roots = rootnodes.backend(arcs, nodes)
  leafs = leafnodes.backend(arcs, nodes)
  to.be.ignored = (arcs[, "from"] %in% roots) | (arcs[, "to"] %in% leafs)

  # compute the loop counter of each (relevant) arc.
  loops = loop.counter(arcs[!to.be.ignored, , drop = FALSE], nodes)

  if (debug) {

    cat("* checking whether the (partially directed) graph is acyclic.\n")

    if (any(to.be.ignored)) {

      cat("  > ignored arcs:", length(which(to.be.ignored)), "\n")
      print(arcs[to.be.ignored, , drop = FALSE])

    }#THEN

    cat("  > loops:\n")
    print(loops)

  }#THEN

  return(sum(loops[, "loops"]) == 0)

}#IS.PARTIALLY.ACYCLIC

# compute the data / cells ratio.
obs.per.cell = function(x, y, z = NULL, data) {

  # return +Inf  for continuous data to bypass countermeasures
  # thought for scarce discrete data.
  if (is.data.continuous(data))
    return(Inf)

  if (is.null(z) || (length(z) == 0)) {

    nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]))

  }#THEN
  else if (is.character(z)) {

    if (length(z) == 1)
      nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) * nlevels(data[,z]))
    else if (length(z) > 1)
      nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) *
        prod(sapply(z, function(col) { nlevels(data[, col]) } )))

  }#THEN
  else if (is.factor(z)) {

    nrow(data) / (nlevels(data[,x]) * nlevels(data[,y]) * nlevels(z))

  }#ELSE

}#OBS.PER.CELL

