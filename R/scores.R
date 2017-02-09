
# compute individual node contributions to the network score.
per.node.score = function(network, data, score, targets, extra.args,
    debug = FALSE) {

  .Call("per_node_score",
        network = network,
        data = data,
        score = score,
        targets = targets,
        extra.args = extra.args,
        debug = debug)

}#PER.NODE.SCORE

# complete a prior over arcs as per Castelo and Siebes.
cs.completed.prior = function(beta, nodes, learning = FALSE) {

  beta = .Call("castelo_completion",
               prior = beta,
               nodes = nodes,
               learning = learning)

  class(beta) = c("prior", "prior.cs", "data.frame")
  attr(beta, "nodes") = nodes

  return(beta)

}#CS.COMPLETED.PRIOR

# compute the optimal imaginary sample size for a discrete network.
alpha.star.backend = function(x, data, debug = FALSE) {

  .Call("alpha_star",
    x = x,
    data = data,
    debug = debug)

}#ALPHA.STAR.BACKEND

