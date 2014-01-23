
# convert a set of arcs to a (real) edge list.
arcs2elist = function(arcs, nodes, weights = NULL, nid = TRUE, sublist = TRUE,
    parents = FALSE) {

  .Call("arcs2elist",
        arcs = arcs,
        nodes = nodes,
        weigths = weights,
        nid = nid,
        sublist = sublist,
        parents = parents)

}#ARCS2ELIST

# convert an edge list into an arc set.
elist2arcs = function(elist) {

  .Call("elist2arcs",
        elist = elist)

}#ELIST2ARCS

