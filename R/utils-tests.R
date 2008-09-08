# check a matrix for symmetry.
is.symmetric = function(m) {

  # kill all the names; identical may return false otherwise.
  colnames(m) = rownames(m) = NULL

  identical(m, t(m))

}#IS.SYMMETRIC

# check whether a graph is completely directed (has no undirected arcs).
is.dag = function(arcs, nodes) {

  if (nrow(arcs) == 0)
    return(TRUE)

  .Call("is_dag",
        arcs = factor(arcs),
        nnodes = length(nodes),
        PACKAGE = "bnlearn")

}#IS.DAG

# alias of is.dag for readable code.
is.pdag = function(arcs, nodes) !is.dag(arcs, nodes)

# generic is.acyclic backend (calls the right one for the graph at hand).
is.acyclic = function(arcs, nodes, debug = FALSE) {

  is.acyclic.backend(arcs = arcs, nodes = nodes,
    directed = is.dag(arcs, nodes), debug = debug)

}#IS.ACYCLIC

# check whether a (partially directed) graph is acyclic.
# Extension to mixed graphs of the proof for proposition 1.4.3
# in "Digraphs Theory, Algorithms and Applications" by
# Bang-Jensen and Gutin, page 13.
is.acyclic.backend = function(arcs, nodes, directed = FALSE, return.nodes = FALSE,
  debug = FALSE) {

  # no arcs, no cycles.
  if (nrow(arcs) == 0)
    if (return.nodes)
      return(character(0))
    else
      return(TRUE)

  if (debug)
    cat("* building the adjacency matrix.\n")

  amat = arcs2amat(arcs, nodes)

  if (debug)
    cat("* checking whether the (directed) graph is acyclic.\n")

  in.loop = nodes

  # to be part of a loop, a node needs to have both a parent
  # and a children; remove (iteratively) all the nodes which
  # do not fit this description.
  for (i in seq(1, length(nodes))) {

    # save a copy of the names of the 'good' nodes at the beginning of each
    # iteration; it will be used to check if any node is removed and avoid
    # useless iterations for cyclic graphs.
    last.known.in.loop = in.loop

    if (debug)
      cat("  > nodes which may be part of a cycle: '", in.loop, "'\n")

    row.sums = rowSums(amat[in.loop, in.loop, drop = FALSE])
    col.sums = colSums(amat[in.loop, in.loop, drop = FALSE])

    # iteratively remove nodes which have either no parents
    # or no children.
    in.loop = intersect(names(row.sums[row.sums > 0]),
                     names(col.sums[col.sums > 0]))

    if (debug)
      cat("  > nodes which may be part of a cycle: '", in.loop, "'\n")

    # the following code is useless for directed graphs; skip.
    if (!directed) {

      # a single undirected arc cannot cause a node to be part of a loop;
      # thic condition can be identified as follows:
      # 1) row total equal to 1 (only one outward arc).
      row.sums = rowSums(amat[in.loop, in.loop, drop = FALSE])
      # 2) column totla equal to 1 (only one inward arc).
      col.sums = colSums(amat[in.loop, in.loop, drop = FALSE])
      # 3) row %*% col = 1 which means that the inward and outward arcs are
      #      one and the same.
      single.undirected = rowSums(amat[in.loop, in.loop] * t(amat[in.loop, in.loop]))

      in.loop = in.loop[!((single.undirected == 1) & (row.sums == 1) & (col.sums == 1))]

    }#THEN

    # a loop requires at least three distinct nodes.
    if (length(in.loop) < 3) {

      if (debug) cat("  @ no cycle found.\n")

      if (return.nodes)
        return(character(0))
      else
        return(TRUE)

    }#THEN

    # if the 'good' nodes are the same as in the last iteration, the graph is
    # cyclic; break the loop and return accordingly.
    if (length(in.loop) == length(last.known.in.loop)) break

  }#FOR

  if (debug)
    cat("  @ these nodes are part of one or more cycles: '", in.loop, "'\n")

  if (return.nodes)
    return(in.loop)
  else
    return(FALSE)

}#IS.ACYCLIC.BACKEND

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

