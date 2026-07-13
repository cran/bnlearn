
# sanity check adjacency matrices.
check.amat = function(amat, nodes) {

  # a node is needed.
  if (missing(amat))
    stop("no adjacency matrix specified.")
  # the adjacency matrix must, well, be a matrix.
  if (!is(amat, "matrix") || (ncol(amat) != nrow(amat)) || (length(dim(amat)) != 2))
    stop("an adjacency matrix must be a 2-dimensional square matrix.")
  # check the dimensions against the number of nodes in the graph.
  if (any(dim(amat) != length(nodes)))
    stop("the dimensions of the adjacency matrix do not agree with the number of nodes in the graph.")
  # column names must be valid node labels.
  if (!is.null(colnames(amat)))
    if (any(colnames(amat) %!in% nodes))
      stop("node (column label) not present in the graph.")
  # column names must be valid node labels.
  if (!is.null(rownames(amat)))
    if (any(rownames(amat) %!in% nodes))
      stop("node (row label) not present in the graph.")
  # column names must match with row names.
  if (!is.null(colnames(amat)) && !is.null(rownames(amat))) {

    if (!identical(colnames(amat), rownames(amat)))
      stop("row/column names mismatch in the adjacency matrix.")

    if (!identical(colnames(amat), nodes) || !identical(rownames(amat), nodes)) {

      warning("rearranging the rows/columns of the adjacency matrix.")

      amat = amat[nodes, nodes, drop = FALSE]

    }#THEN

  }#THEN
  # make really sure the adjacency matrix is made up of integers.
  if (storage.mode(amat) != "integer")
    storage.mode(amat) = "integer"
  # check the elements of the matrix.
  if (!all((amat == 0L) | (amat == 1L)))
    stop("all the elements of an adjacency matrix must be equal to either 0 or 1.")
  # no arcs from a node to itself.
  if (any(diag(amat) != 0))
    stop("the elements on the diagonal must be zero.")

  return(amat)

}#CHECK.AMAT

# sanity check arc lists.
check.alist = function(alist, nodes) {

  # the adjacency list must be a list, with one element for each node.
  if (!is.list(alist))
    stop("the arc list must be a list.")
  if (length(alist) != length(nodes))
    stop("the arc list must have one element for each node in the network.")
  # if the elements are named, the names must match the node labels.
  if (is.null(names(alist))) {

    names(alist) = nodes

  }#THEN
  else {

    if (!setequal(names(alist), nodes))
      stop("the arc list element names and the network node labels do not match.")

    # reorder the elements of the arc list to match the nodes.
    alist = alist[nodes]

  }#ELSE

  # replace NULLS with character(0), to allow shorthand notation.
  are.null = sapply(alist, is.null)
  if (any(are.null))
    alist[are.null] = list(character(0))

  # the elements of the arc list must be character vectors.
  are.strings = sapply(alist, is.string.vector)
  if (any(!are.strings))
    stop("the following elements of the arc list are not character vectors: ",
      paste(names(which(!are.strings)), collapse = " "), ".")
  # the character vectors can only contain node labels.
  are.labels = sapply(alist, function(x) all(x %in% nodes))
  if (any(!are.labels))
    stop("the following elements of the arc list contain invalid node labels: ",
      paste(names(which(!are.labels)), collapse = " "), ".")
  # the elements of the arc list must not contain duplicate strings.
  have.dupes = sapply(alist, anyDuplicated)
  if (any(have.dupes))
    stop("the following elements of the arc list contain duplicated labels: ",
      paste(names(which(have.dupes)), collapse = " "), ".")
  # no arcs from a node to itself.
  loops = sapply(nodes, function(x) x %in% alist[[x]])
  if (any(loops))
    stop("the following nodes have loops:",
      paste(names(which(loops)), collapse = " "), ".")

  return(alist)

}#CHECK.ALIST

