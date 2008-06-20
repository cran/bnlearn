# do a partial ordering of the nodes of a graph.
schedule = function(x, debug = FALSE) {

  nodes = rootnodes(x)
  to.do = rep(0, length(x$nodes))
  names(to.do) = names(x$nodes)

  # this is a very simple implementation of a non-recursive breadth
  # first search; infinite loops are impossible because the longest
  # path in a graph is {#nodes}-long.
  for (step in 1:length(x$nodes)) {

    to.do[nodes] = step

    if (debug) {

      cat("* at depth", step, "the ordering of the nodes is:\n")
      print(to.do)

    }#THEN

    nodes = sapply(nodes, function(node) { children(x, node) }, 
              simplify = FALSE)

    if (is.list(nodes)) 
      nodes = do.call("c", nodes)

    nodes = unique(nodes)

    if (length(nodes) == 0) break

  }#FOR

  return(names(sort(to.do)))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, cap 3.6.1.
rbn.discrete = function(x, n, data, debug = FALSE) {

  to.do = schedule(x, debug = debug)
  result = list()

  for(node in to.do){

    node.parents = parents.backend(x$arcs, node)
    node.levels = levels(data[, node])

    if (debug) {

      cat("* simulating node", node, "with parents '", node.parents, "'.\n")

    }#THEN

    if (length(node.parents) == 0) {

      prob = table(data[, node])/length(data[, node])
      result[[node]] = factor(x = sample(node.levels, n,
          replace = TRUE, prob = prob), levels = node.levels)

    }#THEN
    else {

      # build a table of the original data to compute the conditional
      # probabilities for each configuration of the parents.
      # if there is only a single parent, use it as it is.
      if (length(node.parents) == 1)
        config = data[, node.parents]
      else
        config = configurations(data[,node.parents])

      tab = table(data = data[, node], cfg = config)

      # build an array with the configurations of the simulated data.
      # if there is only a single parents, use it as is.
      if (length(node.parents) == 1)
        config2 = result[[node.parents]]
      else
        config2 = factor(apply(as.data.frame(result)[, node.parents], 1,
                      paste, sep = "", collapse = ":"))

      result[[node]] =
      sapply(config2, function(cfg) {

        factor(x = sample(node.levels, 1, replace = TRUE,
          prob = tab[, cfg]/sum(tab[, cfg])), levels = node.levels)

      })

    }#ELSE


  }#FOR

  as.data.frame(result)[, nodes(x)]

}#RBN.DISCRETE
