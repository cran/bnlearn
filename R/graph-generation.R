
# generate a random directed acyclic graph.
# Original code from the pcalg package, released
# under "GPLv2 or later" with copyright "Markus Kalisch, 2006".
random.graph.backend = function (nodes, prob) {

    n = length(nodes)

    # sample outgoing arcs for the first n - 2 nodes.
    # and, yes, this preseves node ordering.
    edL = lapply(seq(length = n - 2), function(i) {

      outgoing = sample(nodes[seq(i + 1, n)], size = rbinom(1, n - i, prob))

      if (length(outgoing) > 1)
        list(nbr = outgoing)
      else
        list(nbr = character(0))

    })

    # add outgoing arcs for the "n - 1"th node.
    if (rbinom(1, 1, prob) == 1)
        edL[[n - 1]] = list(nbr = nodes[n])
    else
        edL[[n - 1]] = list(nbr = character(0))

    # the last node must no have outgoing arcs.
    edL[[n]] = list(nbr = character(0))

    # set the names of the list elements.
    names(edL) = nodes

    res = empty.graph.backend(nodes)
    res$arcs = mb2arcs(edL, nodes)
    res$nodes = cache.structure(nodes, res$arcs)

return(res)

}#RANDOM.GRAPH.BACKEND

# generate an empty 'bn' object given a set of nodes.
empty.graph.backend = function(nodes) {

  arcs = matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
  nodes.structure = structure(lapply(nodes, function(n) {
        list(mb = character(0), nbr = character(0), parents = character(0),
          children = character(0)) }), names = nodes)
  res = structure(list(
    learning = list(
      nodes = nodes.structure,
      arcs = arcs,
      whitelist = NULL,
      blacklist = NULL,
      test = "none",
      ntests = 0,
      algo = "rnd",
      args = list()
    ),
    nodes = nodes.structure,
    arcs = arcs
  ), class = "bn")

  res

}#EMPTY.GRAPH.BACKEND
