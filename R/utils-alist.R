
# convert a set of arcs to a (real) edge list.
arcs2alist = function(arcs, nodes, weights = NULL, nid = TRUE, sublist = TRUE,
    parents = FALSE) {

  .Call(call_arcs2alist,
        arcs = arcs,
        nodes = nodes,
        weights = weights,
        nid = nid,
        sublist = sublist,
        parents = parents)

}#ARCS2ALIST

# convert an edge list into an arc set.
alist2arcs = function(alist, parents = FALSE) {

  .Call(call_alist2arcs,
        alist = alist,
        parents = parents)

}#ALIST2ARCS

