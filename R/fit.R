
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, method = "mle", extra.args, debug = FALSE) {

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
  fitted = sapply(names(x$nodes), fit, x = x, data = data, method = method,
             extra.args = extra.args, debug = debug, simplify = FALSE)
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

bn.fit.backend.discrete = function(x, node, data, method, extra.args,
    debug = FALSE) {

  # store the labels of the parents and the children to get them only once.
  parents = x$nodes[[node]]$parents
  children = x$nodes[[node]]$children

  if (debug) {

    cat("* fitting parameters of node", node, "(discrete).\n")

    if (length(parents) > 0)
      cat("  > found parents:", parents, "\n")

  }#THEN

  if (method == "mle") {

    # the parameters of the multinomial distribution are the probabilities
    # of the levels of the node and the configurations of its parents.
    tab = table(data[, c(node, parents), drop = FALSE])

  }#THEN
  else {

    # the parameters of the multinomial distribution are the expected values
    # of the corresponding parameters of the posterior Dirichlet distribution.
    tab = table(data[, c(node, parents), drop = FALSE])
    tab = tab + extra.args$iss / prod(dim(tab))

  }#ELSE

  # switch from joint to conditional probabilities.
  tab = prop.table(tab, margin = seq(length(parents) + 1)[-1])
  # this is to preserve the ordering of the factor.
  class = ifelse(is(data[, node], "ordered"), "bn.fit.onode", "bn.fit.dnode")

  # marginal tables have no dimension names in bnlearn.
  tab = cptattr(tab)

  structure(list(node = node, parents = parents, children = children,
    prob = tab), class = class)

}#BN.FIT.BACKEND.DISCRETE

bn.fit.backend.continuous = function(x, node, data, method, extra.args,
    debug = FALSE) {

  # cache the sample size.
  n = nrow(data)
  # store the labels of the parents and the children to get them only once.
  parents = x$nodes[[node]]$parents
  children = x$nodes[[node]]$children

  if (debug)
    cat("* fitting parameters of node", node, "(continuous).\n")

  if (length(parents) == 0) {

    # the only parameters are the intercept of the model (which is the only
    # coefficient) and the standard deviation of the residuals.
    mean = mean(data[, node])
    coefs = c("(Intercept)" = mean)
    resid = data[, node] - mean
    sd = cgsd(data[, node])

    structure(list(node = node, parents = parents, children = children,
      coefficients = coefs, residuals = resid,
      fitted.values = rep(mean, n), sd = sd), class = "bn.fit.gnode")

  }#THEN
  else {

    if (debug)
      cat("  > found parents:", parents, "\n")

    # the relevant quantities of the normal distribution are those usually found
    # in an lm object: coefficients, fitted values and residuals, plus the
    # residuals standard deviation.
    qr.x = qr(minimal.qr.matrix(data, parents))
    coefs = structure(qr.coef(qr.x, data[, node]), names = c("(Intercept)", parents))
    fitted = qr.fitted(qr.x, data[, node])
    resid = qr.resid(qr.x, data[, node])
    sd = cgsd(resid, p = length(parents) + 1)

    structure(list(node = node, parents = parents, children = children,
      coefficients = c(coefs), residuals = resid,
      fitted.values = fitted, sd = sd), class = "bn.fit.gnode")

  }#ELSE

}#BN.FIT.BACKEND.CONTINUOUS

bn.fit.backend.mixedcg = function(x, node, data, method, extra.args,
    debug = FALSE) {

  # cache the node's data.
  node.data = data[, node]

  # discrete nodes are parameterized by CPTs.
  if (is(node.data, "factor"))
    return(bn.fit.backend.discrete(x = x, node = node, data = data,
      method = method, extra.args = extra.args, debug = debug))

  # store the labels of the parents and the children to get them only once.
  parents = x$nodes[[node]]$parents
  children = x$nodes[[node]]$children

  if (length(parents) == 0) {

    # this node has no parents, so the local distribution is just a linear
    # regression against the sample mean.
    return(bn.fit.backend.continuous(x = x, node = node, data = data,
      method = method, extra.args = extra.args, debug = debug))

  }#THEN
  else {

    parents.data = data[, parents, drop = FALSE]
    node.type = data.type(parents.data)

    # if all parents are continuous, the local distribution is just a (single)
    # linear regression with parents as regressors.
    if (node.type == "continuous")
      return(bn.fit.backend.continuous(x = x, node = node, data = data,
        method = method, extra.args = extra.args, debug = debug))
    else {

      continuous.parents = names(which(sapply(parents.data, is.numeric)))
      discrete.parents = parents[parents %!in% continuous.parents]
      configs = configurations(parents.data[discrete.parents])

      if (debug) {

        cat("* fitting parameters of node", node, "(conditional Gaussian).\n")
        cat("  > found parents:", parents, "\n")

      }#THEN

      # perform a linear regression for each configuration of the discrete
      # parents.
      fitted = by(data = data[, c(node, continuous.parents), drop = FALSE],
                   INDICES = configs, FUN = function(x) {

        qr.x = qr(minimal.qr.matrix(x, continuous.parents))
        coefs = structure(qr.coef(qr.x, x[, 1]), names = c("(Intercept)", continuous.parents))
        fitted = qr.fitted(qr.x, x[, 1])
        resid = qr.resid(qr.x, x[, 1])
        sd = cgsd(resid, p = ncol(x))

        # zero all NA coefficients; they arise when a continuous parent is
        # constant for one or more configurations of the discrete parents,
        # but not all of them.
        coefs[is.na(coefs)] = 0

        return(list(coefs = coefs, sd = sd, fitted = fitted, resid = resid))

      })

      residuals = numeric(nrow(data))
      for (cfg in levels(configs))
        residuals[configs == cfg] = fitted[[cfg]]$resid
      fitted.values = numeric(nrow(data))
      for (cfg in levels(configs))
        fitted.values[configs == cfg] = fitted[[cfg]]$fitted

      get.coefficients = function(x) {

        coef.matrix = sapply(x, function(node) {

          if (is.null(node))
            return(structure(rep(NaN, length(continuous.parents) + 1),
              names = c("(Intercept)", continuous.parents)))
          else
            return(node$coefs)

        })

        # make sure the coefficients matrix really is a matrix, sapply may
        # simplify it into a vector.
        if (!is.matrix(coef.matrix))
          coef.matrix = matrix(coef.matrix, nrow = 1, dimnames = list("(Intercept)",
                          as.character(1:length(coef.matrix) - 1)))

        return(coef.matrix)

      }#GET.COEFFICIENTS

      get.sd = function(x) {

        ifelse(is.null(x), NaN, x$sd)

      }#GET.SD

      structure(list(node = node, parents = parents, children = children,
        dparents = which(parents %in% discrete.parents),
        gparents = which(parents %in% continuous.parents),
        dlevels = lapply(parents.data[parents %in% discrete.parents], levels),
        coefficients = get.coefficients(fitted),
        residuals = residuals, fitted.values = fitted.values,
        configs = configs, sd = sapply(fitted, get.sd)
      ), class = "bn.fit.cgnode")

    }#ELSE

  }#ELSE

}#BN.FIT.BACKEND.MIXEDCG

