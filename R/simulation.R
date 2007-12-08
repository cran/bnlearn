schedule = function(x, debug = FALSE) {

to.do = c()

  schedule.node = function(node, x, to.do) {

    if (debug) cat("* checking node", node, "\n")

    node.parents = x$nodes[[node]]$parents
    if (length(node.parents) >= 1) {

      # if a node has parents, schedule them before the node itself.
      if (debug) cat("  > node has parents '", node.parents, "'\n")
      for(parent in node.parents)
        to.do = schedule.node(parent, x, to.do)

      # then schdule the node if it's not already in queue.
      if (!(node %in% to.do))
        to.do = c(to.do, node)
      if (debug) cat("  @ current node ordering is: ", to.do, "\n")

    }#THEN
    else {

      if (!(node %in% to.do)) {

        # if a node has no parents, schedule it.
        to.do = c(node, to.do)
        if (debug) cat("  @ current node ordering is: ", to.do, "\n")

      }#THEN
      else if (debug) {

        cat("  > node", node, "already scheduled.\n")

      }#ELSE

    }#ELSE

    to.do

  }#SCHEDULE.NODE

  for (node in nodes(x))
    to.do = schedule.node(node, x, to.do)

  return(to.do)

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
