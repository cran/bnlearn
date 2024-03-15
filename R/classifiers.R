
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
tan.backend = function(data, training, explanatory, whitelist, blacklist, mi,
    root, debug = FALSE) {

  # cache the node set.
  nodes = c(training, explanatory)
  # create the empty graph.
  res = empty.graph(nodes)
  # create the set of arcs outgoing from the training variable.
  class.arcs = matrix(c(rep(training, length(explanatory)), explanatory),
                 ncol = 2, byrow = FALSE)

  # separate features and target class variable in data and metadata.
  features.data =
    .data.frame.column(data, explanatory, drop = FALSE, keep.names = TRUE)
  features.data = .data.frame(features.data)
  attr(features.data, "metadata") = collect.metadata(features.data)
  class.data = .data.frame.column(data, training, drop = TRUE)

  # call chow-liu to build the rest of the network.
  chow.liu.arcs =
    chow.liu.backend(x = features.data, nodes = explanatory, estimator = mi,
      whitelist = whitelist, blacklist = blacklist, conditional = class.data,
      debug = debug)

  # set the directions of the arcs in the Chow-Liu tree.
  chow.liu.arcs = .Call(call_tree_directions,
                        arcs = chow.liu.arcs,
                        nodes = explanatory,
                        root = root,
                        debug = FALSE)

  # merge learned and predetermined arcs.
  res$arcs = arcs.rbind(class.arcs, chow.liu.arcs)
  # update the network structure.
  res$nodes = cache.structure(nodes, arcs = res$arcs)
  # set a second class "bn.tan" to reroute method dispatch as needed.
  class(res) = c("bn.tan", "bn")

  return(res)

}#TAN.BACKEND

