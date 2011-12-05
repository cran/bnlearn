# backend of the naive bayes classifier.
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
  # set a second class "bn.naive" to reroute the dispatch as needed.
  class(res) = c("bn.naive", "bn")

  return(res)

}#NAIVE.BAYES.BACKEND
