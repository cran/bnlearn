
# reconstruct the equivalence class of a network.
cpdag.backend = function(x, moral = TRUE, fix = FALSE, wlbl = FALSE, debug = FALSE) {

  nodes = names(x$nodes)

  # reset wlbl if x contains no whitelist and no blacklist.
  if ((is.null(x$learning$whitelist)) && (is.null(x$learning$blacklist)))
    wlbl = FALSE

  amat = .Call("cpdag",
               arcs = x$arcs,
               nodes = nodes,
               moral = moral,
               fix = fix || wlbl,
               whitelist = x$learning$whitelist,
               blacklist = x$learning$blacklist,
               illegal = x$learning$illegal,
               debug = debug)

  # update the arcs of the network.
  x$arcs = amat2arcs(amat, nodes)
  # update the network structure.
  x$nodes = cache.structure(nodes, amat = amat)

  return(x)

}#CPDAG.BACKEND

# reconstruct the arc set of the equivalence class of a network.
cpdag.arc.backend = function(nodes, arcs, moral = TRUE, fix.directed = FALSE,
    debug = FALSE) {

  amat = .Call("cpdag",
               arcs = arcs,
               nodes = nodes,
               moral = moral,
               fix = fix.directed,
               whitelist = NULL,
               blacklist = NULL,
               illegal = NULL,
               debug = debug)

  return(amat2arcs(amat, nodes))

}#CPDAG.ARCS.BACKEND

# backend to get a DAG out of a CPDAG (still in the same equivalence class).
cpdag.extension = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # update the arcs of the network.
  x$arcs = cpdag.arc.extension(arcs = x$arcs, nodes = nodes, debug = debug)
  # update the network structure.
  x$nodes = cache.structure(nodes, arcs = x$arcs, debug = debug)

  return(x)

}#CPDAG.EXTENSION

# backend to get a set of directed arcs out of a CPDAG.
cpdag.arc.extension = function(arcs, nodes, debug = FALSE) {

  .Call("pdag_extension",
        arcs = arcs,
        nodes = nodes,
        debug = debug)

}#CPDAG.ARC.EXTENSION

