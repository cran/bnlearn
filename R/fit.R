
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, cluster = NULL, method, extra.args,
    data.info, keep.fitted = TRUE, debug = FALSE) {

  # define the fitting functions.
  if (method %in% c("mle", "bayes"))
    fit = bn.fit.backend.discrete
  else if (method == "mle-g")
    fit = bn.fit.backend.continuous
  else if (method == "mle-cg")
    fit = bn.fit.backend.mixedcg
  else if (method == "hdir")
    fit = bn.fit.backend.hierarchical

  # fit the parameters of each node.
  fitted = smartSapply(cluster, names(x$nodes), fit, dag = x, data = data,
             method = method, extra.args = extra.args, data.info = data.info,
             keep.fitted = keep.fitted, debug = debug)
  # preserve any additional class of the original bn object.
  orig.class = class(x)
  class = c(orig.class[orig.class != "bn"], "bn.fit",
            fitted.from.data[[data.info$type]])
  # preserve the training node label from Bayesian network classifiers.
  if (x$learning$algo %in% classification.algorithms)
    fitted = structure(fitted, class = class, training = x$learning$args$training)
  else
    fitted = structure(fitted, class = class)

  return(fitted)

}#BN.FIT.BACKEND

bn.fit.backend.discrete = function(dag, node, data, method, extra.args,
    data.info, keep.fitted = TRUE, debug = FALSE) {

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

  # the "method" argument is not used explicity, but if "mle" the imaginary
  # sample size is not defined and thus is read as NULL.
  cptable =
    .Call("call_classic_discrete_parameters",
          data = data,
          node = node,
          parents = parents,
          iss = extra.args$iss,
          replace.unidentifiable = isTRUE(extra.args$replace.unidentifiable),
          missing = !all(data.info$complete.nodes[c(node, parents)]))

  # this is to preserve the ordering of the factor.
  class = ifelse(ordered.factor, "bn.fit.onode", "bn.fit.dnode")

  # marginal tables have no dimension names in bnlearn.
  cptable = cptattr(cptable)

  if (debug)
    cat("  > fitted ", length(cptable),
      ifelse(length(parents) > 0, " conditional", " marginal"),
      " probabilities.\n", sep = "")

  structure(list(node = node, parents = parents, children = children,
    prob = cptable), class = class)

}#BN.FIT.BACKEND.DISCRETE

bn.fit.backend.hierarchical = function(dag, node, data, method, extra.args,
    data.info, keep.fitted = TRUE, debug = FALSE) {

  group = extra.args$group
  parents = setdiff(dag$nodes[[node]]$parents, group)
  children = dag$nodes[[node]]$children
  # is the data an ordered or unordered factor?
  ordered.factor = is(data[, node], "ordered")

  if (debug) {

    cat("* fitting parameters of node", node,
      ifelse(ordered.factor, "(ordinal).\n", "(discrete).\n"))

  }#THEN

  # the grouping variable is different from the rest since it has no pooling.
  if (node == group) {

    # fit the parameters of the grouping node using maximum likelihood.
    return(bn.fit.backend.discrete(dag = dag, node = node, data = data,
             data.info = data.info, method = "mle",
             extra.args = list(replace.unidentifiable = FALSE),
             debug = FALSE))

  }#THEN

  if (debug) {

    cat("  > grouping variable:", group, "\n")
    if (length(parents) > 0)
      cat("  > found parents:", parents, "\n")

  }#THEN

  cptable =
    .Call("call_hierarchical_dirichlet_parameters",
          data = data,
          node = node,
          parents = parents,
          group = group,
          alpha0 = extra.args$alpha0,
          iss = extra.args$iss,
          missing = !all(data.info$complete.nodes[c(node, parents, group)]),
          debug = debug)

  # this is to preserve the ordering of the factor.
  class = ifelse(ordered.factor, "bn.fit.onode", "bn.fit.dnode")

  # marginal tables have no dimension names in bnlearn.
  cptable = cptattr(cptable)

  if (debug)
    cat("  > fitted ", length(cptable),
      ifelse(length(parents) > 0, " conditional", " marginal"),
      " probabilities.\n", sep = "")

  structure(list(node = node, parents = c(parents, group), children = children,
                 prob = cptable), class = class)

}#BN.FIT.BACKEND.HIERARCHICAL

bn.fit.backend.continuous = function(dag, node, data, method, extra.args,
    data.info, keep.fitted = TRUE, debug = FALSE) {

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

  # the key elements found in an lm object: coefficients, fitted values and
  # residuals, plus the standard error.
  lmreg =
    .Call("call_gaussian_ols_parameters",
          data = data,
          node = node,
          parents = parents,
          keep.fitted = keep.fitted,
          replace.unidentifiable = isTRUE(extra.args$replace.unidentifiable),
          missing = !all(data.info$complete.nodes[c(node, parents)]))

  structure(c(list(node = node, parents = parents, children = children),
    lmreg), class = "bn.fit.gnode")

}#BN.FIT.BACKEND.CONTINUOUS

bn.fit.backend.mixedcg = function(dag, node, data, method, extra.args,
    data.info, keep.fitted = TRUE, debug = FALSE) {

  # cache the node's data.
  node.data = data[, node]

  # discrete nodes are parameterized by CPTs.
  if (is(node.data, "factor"))
    return(bn.fit.backend.discrete(dag = dag, node = node, data = data,
      method = "mle", extra.args = extra.args, data.info = data.info,
      debug = debug))

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children

  if (length(parents) == 0) {

    # this node has no parents, so the local distribution is just a linear
    # regression against the sample mean.
    return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
      method = "mle-g", extra.args = extra.args, keep.fitted = keep.fitted,
      data.info = data.info, debug = debug))

  }#THEN
  else {

    parents.data = data[, parents, drop = FALSE]
    type.from.parents = data.type(parents.data)

    # if all parents are continuous, the local distribution is just a (single)
    # linear regression with parents as regressors.
    if (type.from.parents == "continuous") {

      return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
        method = "mle-g", extra.args = extra.args, keep.fitted = keep.fitted,
        data.info = data.info, debug = debug))

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
      lmreg =
        .Call("call_mixture_gaussian_ols_parameters",
              data = data,
              node = node,
              parents = continuous.parents,
              configs = configs,
              keep.fitted = keep.fitted,
              replace.unidentifiable = isTRUE(extra.args$replace.unidentifiable),
              missing = !all(data.info$complete.nodes[c(node, parents)]))

      structure(c(list(node = node, parents = parents, children = children,
        dparents = which(parents %in% discrete.parents),
        gparents = which(parents %in% continuous.parents),
        dlevels = lapply(parents.data[parents %in% discrete.parents], levels)),
        lmreg), class = "bn.fit.cgnode")

    }#ELSE

  }#ELSE

}#BN.FIT.BACKEND.MIXEDCG

