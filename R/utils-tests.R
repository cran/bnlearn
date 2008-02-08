# check a matrix for symmetry.
is.symmetric = function(m) {

  # kill all the names; identical may return false otherwise.
  colnames(m) = rownames(m) = NULL

  identical(m, t(m))

}#IS.SYMMETRIC

# check whether a graph is acyclic.
is.acyclic = function(arcs, nodes) {

  # no arcs, no cycles.
  if (nrow(arcs) == 0) return(TRUE)

  amat = arcs2amat(arcs, nodes)

  all(apply(arcs, 1, function(arc) {

    !has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)

  }))

}#IS.ACYCLIC

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

