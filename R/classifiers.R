
# backend of the naive Bayes classifier.
naive.bayes.backend = function(data, training, explanatory) {

  # cache the node set.
  nodes = c(training, explanatory)
  # create the empty graph.
  res = empty.graph(nodes)
  # create the set of arcs outgoing from the training variable.
  res$arcs = matrix(c(rep(training, length(explanatory)), explanatory),
               ncol = 2, byrow = FALSE)
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = res$arcs)
  # set a second class "bn.naive" to reroute method dispatch as needed.
  class(res) = c("bn.naive", "bn")

  return(res)

}#NAIVE.BAYES.BACKEND

# backend of the TAN algorithm.
tan.backend = function(data, training, explanatory, whitelist, blacklist,
    mi, root, debug) {

  # set a dummy estimator variable.
  estimator = 1L
  # cache the node set.
  nodes = c(training, explanatory)
  # create the empty graph.
  res = empty.graph(nodes)
  # create the set of arcs outgoing from the training variable.
  class.arcs = matrix(c(rep(training, length(explanatory)), explanatory),
               ncol = 2, byrow = FALSE)

  # call chow-liu to build the rest of the network.
  chow.liu.arcs = chow.liu.backend(x = minimal.data.frame.column(data, explanatory),
                    nodes = explanatory, estimator = estimator,
                    whitelist = whitelist, blacklist = blacklist,
                    conditional = minimal.data.frame.column(data, training, drop = TRUE),
                    debug = debug)
  # set the directions of the arcs in the Chow-Liu tree.
  chow.liu.arcs = .Call("tree_directions",
                        arcs = chow.liu.arcs,
                        nodes = explanatory,
                        root = root,
                        debug = FALSE,
                        PACKAGE = "bnlearn")

  # merge learned and predetermined arcs.
  res$arcs = arcs.rbind(class.arcs, chow.liu.arcs)
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = res$arcs)
  # set a second class "bn.tan" to reroute method dispatch as needed.
  class(res) = c("bn.tan", "bn")

  return(res)

}#TAN.BACKEND

