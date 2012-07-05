# do a partial ordering of the nodes of a graph.
schedule = function(x, start = NULL, reverse = FALSE, debug = FALSE) {

  if (is.null(start))
    nodes = root.leaf.nodes(x, leaf = reverse)
  else
    nodes = start

  to.do = .Call("schedule",
                bn = x,
                root.nodes = nodes,
                reverse = reverse,
                debug = debug,
                PACKAGE = "bnlearn")

  if (is.null(start))
    return(names(sort(to.do)))
  else
    return(names(sort(to.do[to.do > 0])))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, chap 3.6.1.
rbn.discrete = function(x, n, data, debug = FALSE) {

  # fit the bayesian network if needed.
  if (is(x, "bn"))
    fitted = bn.fit.backend(x, data, debug = FALSE)
  else
    fitted = x

  .Call("rbn_discrete",
        fitted = fitted,
        n = as.integer(n),
        debug = debug,
        PACKAGE = "bnlearn")

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
