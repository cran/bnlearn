
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, cluster = NULL, method = "mle", extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # check which type of data we are dealing with.
  type = data.type(data)

  # define the fitting functions.
  if (type %in% discrete.data.types)
    fit = bn.fit.backend.discrete
  else if (type == "continuous")
    fit = bn.fit.backend.continuous
  else if (type == "mixed-cg")
    fit = bn.fit.backend.mixedcg

  # fit the parameters of each node.
  fitted = smartSapply(cluster, names(x$nodes), fit, dag = x, data = data,
             method = method, extra.args = extra.args, keep.fitted = keep.fitted,
             debug = debug)
  # preserve any additional class of the original bn object.
  orig.class = class(x)
  class = c(orig.class[orig.class != "bn"], "bn.fit", fitted.from.data[[type]])
  # preserve the training node label from Bayesian network classifiers.
  if (x$learning$algo %in% classifiers)
    fitted = structure(fitted, class = class, training = x$learning$args$training)
  else
    fitted = structure(fitted, class = class)

  return(fitted)

}#BN.FIT.BACKEND

bn.fit.backend.discrete = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children
  # is the data an ordered or unordered factor?
  ordered.factor = is(data[, node], "ordered")

  if (debug) {

    cat("* fitting parameters of node", node,
      ifelse(ordered.factor, "(ordinal).\n", "(discrete).\n"))

    if (length(parents) > 0)
      cat("  > found parents:", parents, "\n")

  }#THEN

  if (method == "mle") {

    # the parameters of the multinomial distribution are the probabilities
    # of the levels of the node and the configurations of its parents.
    tab = minimal.table(data[, c(node, parents), drop = FALSE])

  }#THEN
  else {

    # the parameters of the multinomial distribution are the expected values
    # of the corresponding parameters of the posterior Dirichlet distribution.
    tab = minimal.table(data[, c(node, parents), drop = FALSE])
    tab = tab + extra.args$iss / prod(dim(tab))

  }#ELSE

  # switch from joint to conditional probabilities.
  tab = prop.table(tab, margin = seq(length(parents) + 1)[-1])
  # this is to preserve the ordering of the factor.
  class = ifelse(ordered.factor, "bn.fit.onode", "bn.fit.dnode")

  # marginal tables have no dimension names in bnlearn.
  tab = cptattr(tab)

  if (debug)
    cat("  > fitted ", length(tab),
      ifelse(length(parents) > 0, " conditional", " marginal"),
      " probabilities.\n", sep = "")

  structure(list(node = node, parents = parents, children = children,
    prob = tab), class = class)

}#BN.FIT.BACKEND.DISCRETE

bn.fit.backend.continuous = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # cache the sample size.
  n = nrow(data)
  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children

  if (debug) {

    cat("* fitting parameters of node", node, "(continuous).\n")
    if (length(parents) > 0)
      cat("  > found parents:", parents, "\n")
    cat("  > fitted", length(parents) + 1, "regression coefficient(s) and",
        "1 standard error.\n")

  }#THEN

  # the relevant quantities of the normal distribution are those usually found
  # in an lm object: coefficients, fitted values and residuals, plus the
  # residuals standard deviation.

  structure(c(list(node = node, parents = parents, children = children),
    fast.lm(data, node, parents, keep.fitted = keep.fitted)),
    class = "bn.fit.gnode")

}#BN.FIT.BACKEND.CONTINUOUS

bn.fit.backend.mixedcg = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # cache the node's data.
  node.data = data[, node]

  # discrete nodes are parameterized by CPTs.
  if (is(node.data, "factor"))
    return(bn.fit.backend.discrete(dag = dag, node = node, data = data,
      method = method, extra.args = extra.args, debug = debug))

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children

  if (length(parents) == 0) {

    # this node has no parents, so the local distribution is just a linear
    # regression against the sample mean.
    return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
      method = method, extra.args = extra.args, keep.fitted = keep.fitted,
      debug = debug))

  }#THEN
  else {

    parents.data = data[, parents, drop = FALSE]
    node.type = data.type(parents.data)

    # if all parents are continuous, the local distribution is just a (single)
    # linear regression with parents as regressors.
    if (node.type == "continuous") {

      return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
        method = method, extra.args = extra.args, keep.fitted = keep.fitted,
        debug = debug))

    }#THEN
    else {

      continuous.parents = names(which(sapply(parents.data, is.numeric)))
      discrete.parents = parents[parents %!in% continuous.parents]
      configs = configurations(parents.data[discrete.parents])

      if (debug) {

        cat("* fitting parameters of node", node, "(conditional Gaussian).\n")
        cat("  > found continuous parents:", continuous.parents, "\n")
        cat("  > found discrete parents:", discrete.parents, "\n")
        cat("  > fitting", nlevels(configs), "x", length(continuous.parents),
            "regression coefficients and", nlevels(configs), "standard errors.\n")

      }#THEN

      # perform a linear regression for each configuration of the discrete
      # parents.
      fitted = fast.cglm(data, node, parents = continuous.parents,
                 configs = configs, keep.fitted = keep.fitted)

      if (keep.fitted) {

        structure(c(list(node = node, parents = parents, children = children,
          dparents = which(parents %in% discrete.parents),
          gparents = which(parents %in% continuous.parents),
          dlevels = lapply(parents.data[parents %in% discrete.parents], levels)),
          fitted), class = "bn.fit.cgnode")

      }#THEN
      else {

        structure(c(list(node = node, parents = parents, children = children,
          dparents = which(parents %in% discrete.parents),
          gparents = which(parents %in% continuous.parents),
          dlevels = lapply(parents.data[parents %in% discrete.parents], levels)),
          fitted), class = "bn.fit.cgnode")

      }#ELSE

    }#ELSE

  }#ELSE

}#BN.FIT.BACKEND.MIXEDCG

