
# convert a set of arcs to a (real) edge list.
arcs2elist = function(arcs, nodes, id = TRUE, sublist = TRUE) {

  .Call("arcs2elist",
        arcs = arcs,
        nodes = nodes,
        id = id,
        sublist = sublist,
        PACKAGE = "bnlearn")

}#ARCS2ELIST

# convert an edge list into an arc set.
elist2arcs = function(elist) {

  .Call("elist2arcs",
        elist = elist,
        PACKAGE = "bnlearn")

}#ELIST2ARCS

