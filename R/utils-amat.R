
# convert a set of arcs to a (real) adjacency matrix.
arcs2amat = function(arcs, nodes) {

  .Call("arcs2amat",
        arcs = as.character(arcs),
        nodes = as.character(nodes),
        PACKAGE = "bnlearn")

}#ARCS2AMAT

# convert an adjacency matrix back to a set of arcs.
amat2arcs = function(a, nodes) {

  .Call("amat2arcs",
        amat = a,
        nodes = nodes,
        PACKAGE = "bnlearn")

}#AMAT2ARCS

