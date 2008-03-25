# convert a set of markov blankets to a (sort of) adjacency matrix.
mb2amat = function(mb, redux = FALSE) {

  if (redux)
    sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]])})
  else
    sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$mb)})

}#MB2AMAT

# convert a set of neighbourhoods to a (real) adjacency matrix.
nbr2amat = function(mb) {

  sapply(names(mb), function(x) {as.numeric(names(mb) %in% mb[[x]]$nbr)})

}#NBR2AMAT

# convert a set of arcs to a (real) adjacency matrix.
arcs2amat = function(arcs, nodes) {

  .Call("arcs2amat",
        arcs = as.character(arcs),
        nodes = as.character(nodes),
        PACKAGE = "bnlearn")

}#ARCS2AMAT

# convert an adjacency matrix back to a set of arcs.
amat2arcs = function(a, nodes, debug = FALSE) {

  .Call("amat2arcs",
        amat = as.integer(a),
        nodes = as.character(nodes),
        PACKAGE = "bnlearn")

}#AMAT2ARCS

