# convert a set of markov blankets to a (sort of) adjacency matrix.
mb2amat = function(mb, stochastic = FALSE, redux = FALSE) {

  if (redux)
    amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]])})
  else
    amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$mb)})

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#MB2AMAT

# convert a set of neighbourhoods to a (real) adjacency matrix.
nbr2amat = function(mb, stochastic = FALSE) {

  amat = sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$nbr)})

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#NBR2AMAT

# convert a set of arcs to a (real) adjacency matrix.
arcs2amat = function(arcs, nodes, stochastic = FALSE) {

  amat = matrix(rep(0, length(nodes)^2), ncol = length(nodes))
  colnames(amat) = rownames(amat) = nodes

  if (nrow(arcs) > 0)
    for (i in 1:nrow(arcs)) {

      amat[arcs[i, "from"], arcs[i, "to"]] = 1

    }#FOR

  if (stochastic)
    amat = t(apply(amat, 1, function(r) { if(sum(r) != 0) r / sum(r) else r}))

  amat

}#ARCS2AMAT

# convert an adjacency matrix back to a set of arcs.
amat2arcs = function(a, nodes, debug = FALSE) {

  # do not panic if there are no arcs.
  if (sum(a) == 0) {

    if (debug) cat("* no arcs in the graph.\n")

    return(matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to"))))

  }#THEN

  if (is.null(colnames(a)))
    colnames(a) = nodes
  if (is.null(rownames(a)))
    rownames(a) = nodes

  arcs = do.call(rbind,
    sapply(rownames(a), function(node) {

      tos = colnames(a)[a[node, ] > 0]

      if (debug) {

        cat("* preocessing node", node, "\n")
        cat("  > adding arcs from", node, "to its parents: '", tos, "'.\n")

      }#THEN

      # something with two colums should be returned anyway, even when
      # there are no actual arcs for that node.
      if (length(tos) > 0)
        cbind(rep(node, length(tos)), tos)
      else
        matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))

    }, simplify = FALSE)
  )

  # add the column names.
  colnames(arcs) = c("from", "to")

  arcs

}#AMAT2ARCS


