# do a partial ordering of the nodes of a graph.
schedule = function(x, debug = FALSE) {

  nodes = root.leaf.nodes(x, leaf = FALSE)

  to.do = .Call("schedule",
                bn = x,
                root.nodes = nodes,
                debug = debug,
                PACKAGE = "bnlearn")

  return(names(sort(to.do)))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, chap 3.6.1.
rbn.discrete = function(x, n, data, debug = FALSE) {

  to.do = schedule(x)
  result = list()

  if (debug)
    cat("* partial node ordering is:", to.do, "\n")

  for (node in to.do) {

    node.parents = x$nodes[[node]]$parents
    node.levels = levels(data[, node])

    if (debug)
      cat("* simulating node", node, "with parents '", node.parents, "'.\n")

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

      # if there are configurations in the generated data which were
      # not observed in the original data; notify the user and print
      # a human readable comparison when in debug mode.
      generated.configurations = as.character(unique(config2))
      observed.configurations = as.character(unique(config))

      if (!all(generated.configurations %in% observed.configurations)) {

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

# a modified Logic Sampling (LS) algorithm for Gaussian data.
rbn.continuous = function(x, n, data, debug = FALSE) {

  to.do = schedule(x)
  result = list()

  intercept = rep(1, nrow(data))

  if (debug)
    cat("* partial node ordering is:", to.do, "\n")

  for (node in to.do) {

    node.parents = x$nodes[[node]]$parents

    if (debug) 
      cat("* simulating node", node, "with parents '", node.parents, "'.\n")

    if (length(node.parents) == 0) {

      mean = mean(data[, node])
      sd = sd(data[, node])

      if (debug)
        cat("  > node", node, "has mean", mean, "and standard deviation", sd, ".\n")

      result[[node]] = rnorm(n, mean, sd)

    }#THEN
    else {

      # compute the regression coefficient and the standard deviation of
      # the residuals with a QR decomposition.
      qr.x = qr(cbind(intercept, data[, node.parents]))
      coefs = qr.coef(qr.x, data[, node])
      rsd = sd(qr.resid(qr.x, data[, node]))

      # compute the predicted values for the data previously generated for
      # the parents of this node.
      names(coefs) = c("intercept", node.parents)
      mean = rep(coefs["intercept"], n)
      for (parent in node.parents)
        mean = mean + result[[parent]] * coefs[parent]

      if (debug) {

        f = paste("(", format(coefs), ")", c("", node.parents),
              sep = "", collapse = " + ")

        cat("  > node", node, "is fitted as\n   ", f, "\n")
        cat("  > residual standard deviation is", rsd, ".\n")

      }

      # add the gaussian noise, and be done.
      result[[node]] = mean + rnorm(n, 0, rsd)

    }#ELSE

  }#FOR

  as.data.frame(result)[, nodes(x)]

}#RBN.CONTINUOUS
