
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, cluster = NULL, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  if (grepl("-em", method)) {

    fitted = bn.fit.backend.hard.em(x = x, data = data, cluster = cluster,
               extra.args = extra.args, keep.fitted = keep.fitted,
               debug = debug)

  }#THEN
  else {

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
               method = method, extra.args = extra.args,
               keep.fitted = keep.fitted, debug = debug)

  }#ELSE

  # preserve any additional class of the original bn object.
  orig.class = class(x)
  class = c(orig.class[orig.class != "bn"], "bn.fit",
            fitted.from.method[[method]])
  # preserve the training node label from Bayesian network classifiers.
  if (x$learning$algo %in% classification.algorithms)
    fitted = structure(fitted, class = class, training = x$learning$args$training)
  else
    fitted = structure(fitted, class = class)

  return(fitted)

}#BN.FIT.BACKEND

# maximum likelihood and posterior parameter estimation for discrete networks.
bn.fit.backend.discrete = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children
  # is the data an ordered or unordered factor?
  ordered.factor = is(data[, node], "ordered")
  # which nodes do not have missing values?
  complete.nodes = attr(data, "metadata")$complete.nodes

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
          missing = !all(complete.nodes[c(node, parents)]))

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

# hierarchical Bayesian parameter estimation for discrete networks.
bn.fit.backend.hierarchical = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  group = extra.args$group
  parents = setdiff(dag$nodes[[node]]$parents, group)
  children = dag$nodes[[node]]$children
  # is the data an ordered or unordered factor?
  ordered.factor = is(data[, node], "ordered")
  # which nodes do not have missing values?
  complete.nodes = attr(data, "metadata")$complete.nodes

  if (debug) {

    cat("* fitting parameters of node", node,
      ifelse(ordered.factor, "(ordinal).\n", "(discrete).\n"))

  }#THEN

  # the grouping variable is different from the rest since it has no pooling.
  if (node == group) {

    # fit the parameters of the grouping node using maximum likelihood.
    return(bn.fit.backend.discrete(dag = dag, node = node, data = data,
             method = "mle", extra.args = list(replace.unidentifiable = FALSE),
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
          missing = !all(complete.nodes[c(node, parents, group)]),
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

# maximum likelihood parameter estimation for Gaussian networks.
bn.fit.backend.continuous = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children
  # which nodes do not have missing values?
  complete.nodes = attr(data, "metadata")$complete.nodes

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
          missing = !all(complete.nodes[c(node, parents)]))

  structure(c(list(node = node, parents = parents, children = children),
    lmreg), class = "bn.fit.gnode")

}#BN.FIT.BACKEND.CONTINUOUS

# maximum likelihood parameter estimation for conditional Gaussian networks.
bn.fit.backend.mixedcg = function(dag, node, data, method, extra.args,
    keep.fitted = TRUE, debug = FALSE) {

  # which nodes do not have missing values?
  complete.nodes = attr(data, "metadata")$complete.nodes
  # cache the node's data.
  node.data = data[, node]

  # discrete nodes are parameterized by CPTs.
  if (is(node.data, "factor"))
    return(bn.fit.backend.discrete(dag = dag, node = node, data = data,
      method = "mle", extra.args = extra.args, debug = debug))

  # store the labels of the parents and the children to get them only once.
  parents = dag$nodes[[node]]$parents
  children = dag$nodes[[node]]$children

  if (length(parents) == 0) {

    # this node has no parents, so the local distribution is just a linear
    # regression against the sample mean.
    return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
      method = "mle-g", extra.args = extra.args, keep.fitted = keep.fitted,
      debug = debug))

  }#THEN
  else {

    parents.data = data[, parents, drop = FALSE]
    type.from.parents = data.type(parents.data)

    # if all parents are continuous, the local distribution is just a (single)
    # linear regression with parents as regressors.
    if (type.from.parents == "continuous") {

      return(bn.fit.backend.continuous(dag = dag, node = node, data = data,
        method = "mle-g", extra.args = extra.args, keep.fitted = keep.fitted,
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
      lmreg =
        .Call("call_mixture_gaussian_ols_parameters",
              data = data,
              node = node,
              parents = continuous.parents,
              configs = configs,
              keep.fitted = keep.fitted,
              replace.unidentifiable = isTRUE(extra.args$replace.unidentifiable),
              missing = !all(complete.nodes[c(node, parents)]))

      structure(c(list(node = node, parents = parents, children = children,
        dparents = which(parents %in% discrete.parents),
        gparents = which(parents %in% continuous.parents),
        dlevels = lapply(parents.data[parents %in% discrete.parents], levels)),
        lmreg), class = "bn.fit.cgnode")

    }#ELSE

  }#ELSE

}#BN.FIT.BACKEND.MIXEDCG

# hard EM implementation following Koller and Friedman.
bn.fit.backend.hard.em = function(x, data, cluster, extra.args, keep.fitted,
    debug) {

  # initial model for imputation.
  if (!is.null(extra.args$start)) {

    current.fitted = extra.args$start

  }#THEN
  else {

    current.fitted =
      bn.fit.backend(x, data = data, cluster = cluster, method = extra.args$fit,
        extra.args = extra.args$fit.args, keep.fitted = FALSE, debug = FALSE)

  }#ELSE

  # set the log-likelihood of the initial model to -Inf, instead of imputing the
  # data and computing it, to make sure to complete at least one iteration.
  current.loglik = -Inf

  # start iterating.
  iter = 1

  while (iter <= extra.args$max.iter) {

    if (debug)
      cat("* expectation-maximization iteration", iter, ".\n")

    # E-step: complete the data by imputing the missing values.
    completed =
      impute.backend(fitted = current.fitted, data = data, cluster = cluster,
        method = extra.args$impute, extra.args = extra.args$impute.args,
        debug = FALSE)

    # check whether the imputation has been successful so that there are no
    # more missing values.
    successfully.imputed = complete.cases(completed)

    if (!all(successfully.imputed))
      stop("impossible to impute observation(s) ",
           paste(which(!successfully.imputed), collapse = ", "), ".")

    # M-step: fit the parameters from the completed data.
    new.fitted =
      bn.fit.backend(x, data = completed, cluster = cluster,
         method = extra.args$fit, extra.args = extra.args$fit.args,
         keep.fitted = keep.fitted, debug = FALSE)

    # compute the node-average log-likelihood with the new parameters, either
    # for the completed validation data set or the completed training ones.
    if (!is.null(extra.args$newdata)) {

      completed =
        impute.backend(fitted = current.fitted, data = extra.args$newdata,
          cluster = cluster, method = extra.args$impute,
          extra.args = extra.args$impute.args, debug = FALSE)

    }#THEN

    new.loglik = loglikelihood(new.fitted, data = completed,
                   propagate.missing = FALSE) / nrow(completed)

    # check convergence: are parameters different enough (if the two network
    # structures are the same, which is not a given in the first iteration)?
    if ((iter > 1) ||
        isTRUE(equal.backend.bn(bn.net(new.fitted), bn.net(current.fitted)))) {

      reldiff.params =
        sapply(names(new.fitted), function(node) {

          diff = abs(.coefficients(new.fitted)[[node]] -
                     .coefficients(current.fitted)[[node]])

          scale = .coefficients(new.fitted)[[node]]
          scale[scale == 0] = sqrt(.Machine$double.eps)

          return(max(diff / scale, na.rm = TRUE))

        })

      reldiff.params = max(reldiff.params)

    }#THEN
    else {

      reldiff.params = Inf

    }#ELSE

    # check convergence: is the new node-average log-likelihood better enough?
    reldiff.loglik = robust.score.difference(new.loglik, current.loglik) /
                       ifelse(is.finite(new.loglik), abs(new.loglik), 1)

    if (debug) {

      cat("  > the relative difference in the parameters is:",
          reldiff.params, ".\n")
      cat("  > the relative difference in log-likelihood is:",
          reldiff.loglik, ".\n")
      if (reldiff.params <= extra.args$params.threshold)
        cat("  @ the difference in the parameters is smaller than the threshold, stopping.\n")
      if (reldiff.loglik <= extra.args$loglik.threshold)
        cat("  @ the difference in log-likelihood is smaller than the threshold, stopping.\n")

    }#THEN

    # stop if either convergence criterion is satisfied.
    if (reldiff.params <= extra.args$params.threshold)
      break
    if (reldiff.loglik <= extra.args$loglik.threshold)
      break

    # update the model and the associated log-likelihood for the next iteration.
    current.fitted = new.fitted
    current.loglik = new.loglik
    # increase the iteration counter.
    iter = iter + 1

  }#WHILE

  if (debug) {

    if ((reldiff.loglik > extra.args$loglik.threshold) &&
        (reldiff.params > extra.args$params.threshold) &&
        (iter == extra.args$max.iter))
      cat("@ reached the maximum number of iterations, stopping.\n")

  }#THEN

  return(current.fitted)

}#BN.FIT.BACKEND.HARD.EM
