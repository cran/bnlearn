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

  # fit the bayesian network if needed.
  if (is(x, "bn"))
    fitted = bn.fit.backend(x, data, debug = FALSE)
  else
    fitted = x

  # schedule the nodes and prepare the return value.
  to.do = schedule(x)
  result = vector(length(fitted), mode = "list")
  names(result) = names(fitted)

  if (debug)
    cat("* partial node ordering is:", to.do, "\n")

  for (node in to.do) {

    node.parents = fitted[[node]]$parents
    node.levels = labels(fitted[[node]]$prob)[[1]]
    prob = fitted[[node]]$prob

    if (debug)
      cat("* simulating node", node, "with parents '", node.parents, "'.\n")

    if (length(node.parents) == 0) {

      result[[node]] = factor(x = sample(node.levels, n,
          replace = TRUE, prob = prob), levels = node.levels)

    }#THEN
    else {

      # get the contingency table of this node against the configurations
      # of its parents.
      config2 = configurations(result[node.parents])
      tab = collapse.table(fitted[[node]]$prob)

      # initialize the vectors to hold the generated data and the subset identifier.
      temp.gen = character(n)
      to.be.generated = logical(n)
      char.configurations = as.character(config2)

      # generate each value according to the right configuration. Iterate on
      # the latter to achieve a reasonable speed.
      for (cfg in unique(char.configurations)) {

        to.be.generated = (char.configurations == cfg)

        if (all(!is.nan(tab[, cfg]))) {

          temp.gen[to.be.generated] = sample(node.levels, length(which(to.be.generated)),
             replace = TRUE, prob = tab[, cfg])

        }#THEN
        else {

          warning(paste("some configurations of the parents of", node,
            "are not present in the original data. NAs will be generated.",
            collpase = "", sep = " "))

          # the conditional probability distribution in this case is unknown;
          # generate a missing value (NA) insted of gthrowing an error.
          temp.gen[to.be.generated] = NA

        }#ELSE

      }#FOR

      # convert the generated values into a factor object.
      result[[node]] = factor(temp.gen, levels = node.levels)

    }#ELSE

  }#FOR

  # WARNING: in-place modification of result to avoid copying potentially
  # huge data sets around.
  minimal.data.frame(result)

  return(result)

}#RBN.DISCRETE

# a modified Logic Sampling (LS) algorithm for Gaussian data.
rbn.continuous = function(x, n, data, debug = FALSE) {

  if (is(x, "bn"))
    fitted = bn.fit.backend(x, data, debug = FALSE)
  else
    fitted = x

  # schedule the nodes and prepare the return value.
  to.do = schedule(x)
  result = vector(length(fitted), mode = "list")
  names(result) = names(fitted)

  if (debug)
    cat("* partial node ordering is:", to.do, "\n")

  for (node in to.do) {

    node.parents = fitted[[node]]$parents

    if (debug)
      cat("* simulating node", node, "with parents '", node.parents, "'.\n")

    if (length(node.parents) == 0) {

      # extract the mean and the standard deviation.
      mean = fitted[[node]]$coefficients
      sd = fitted[[node]]$sd

      if (debug)
        cat("  > node", node, "has mean", mean, "and standard deviation", sd, ".\n")

      result[[node]] = rnorm(n, mean, sd)

    }#THEN
    else {

      # extract the regression coefficients and the standard deviation.
      coefs = fitted[[node]]$coefficients
      rsd = fitted[[node]]$sd

      # compute the predicted values for the data previously generated for
      # the parents of this node.
      mean = rep(coefs["(Intercept)"], n)
      for (parent in node.parents)
        mean = mean + result[[parent]] * coefs[parent]

      if (debug) {

        f = paste("(", format(coefs), ")", c("", node.parents),
              sep = "", collapse = " + ")

        cat("  > node", node, "is fitted as\n   ", f, "\n")
        cat("  > residual standard deviation is", rsd, ".\n")

      }#THEN

      # add the gaussian noise, and be done.
      result[[node]] = mean + rnorm(n, 0, rsd)

    }#ELSE

  }#FOR

  # WARNING: in-place modification of result to avoid copying potentially
  # huge data sets around.
  minimal.data.frame(result)

  return(result)

}#RBN.CONTINUOUS
