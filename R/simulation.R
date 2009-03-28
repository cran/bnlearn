# do a partial ordering of the nodes of a graph.
schedule = function(x, debug = FALSE) {

  nodes = rootnodes.backend(x$arcs, names(x$nodes))
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

    nodes = .Call("schedule_children",
            graph  = x,
            nodes = nodes,
            PACKAGE = "bnlearn")

    nodes = unique(nodes)

    if (length(nodes) == 0) break

  }#FOR

  return(names(sort(to.do)))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, chap 3.6.1.
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

      # build an array with the configurations of the simulated data.
      # if there is only a single parents, use it as is.
      if (length(node.parents) == 1)
        config2 = result[[node.parents]]
      else
        config2 = configurations(as.data.frame(result)[, node.parents])

      # build a table of the original data to compute the conditional
      # probabilities for each configuration of the parents.
      # if there is only a single parent, use it as it is.
      if (length(node.parents) == 1)
        config = data[, node.parents]
      else
        config = configurations(data[, node.parents])

      # add the configurations present in the generated data to the ones
      # observed in the original data; otherwise there might be a column
      # mismatch in the sapply call below.
      levels(config) = union(levels(config), levels(config2))

      # generate the contingency table of the node against the
      # configurations of its parents.
      tab = table(data = data[, node], cfg = config)

      # there are configurations in the generated data which were
      # not observed in the original data; notify the user and print
      # a human readable comparison when in debug mode.
      generated.configurations = as.character(unique(config2))
      observed.configurations = as.character(unique(config))

      if (!setequal(observed.configurations, generated.configurations)) {

        if (debug) {

          # generate a human readable version of the two sets of configurations
          # for the debugging output.
          obc = unique(apply(data[, node.parents], 1, paste, sep = "",
                     collapse = ":"))
          gnc = unique(apply(as.data.frame(result)[, node.parents], 1,
                     paste, sep = "", collapse = ":"))

          cat("  > observed configurations of the parents:\n")
          print(sort(obc))
          cat("  > configurations present in the generated data:\n")
          print(sort(gnc))
          cat("  > configurations not present in the original data:\n")
          print(sort(union( setdiff(obc, gnc), setdiff(gnc, obc))))

        }#THEN

        warning(paste("some configurations of the parents of", node,
          "are not present in the original data. NA's may be generated.",
           collpase = "", sep = " "))

      }#THEN

      # initialize the vector to hold the generated data.
      temp.gen = character(n)
      char.configurations = as.character(config2)

      # generate each value according to its parents' configuration.
      for (cfg in generated.configurations) {

        to.be.generated = (char.configurations == cfg)

        if (sum(tab[, cfg]) != 0) {

          temp.gen[to.be.generated] = sample(node.levels,
             length(which(to.be.generated)), replace = TRUE,
             prob = tab[, cfg]/sum(tab[, cfg]))

        }#THEN
        else {

          # the probability distribution in this case is unknown;
          # generate a missing value (NA).
          temp.gen[to.be.generated] = NA

        }#ELSE

      }#FOR

      # convert the generated values into a factor object.
      result[[node]] = factor(temp.gen, levels = node.levels)

    }#ELSE

  }#FOR

  as.data.frame(result)[, nodes(x)]

}#RBN.DISCRETE
