
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, method = "mle", extra.args, debug = FALSE) {

  n = nrow(data)

  # define the fitting functions.
  if (is.data.discrete(data)) {

    fit = function(node) {

      # store the labels of the parents and the children to get
      # them only once.
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

      # switch from the joint probabilities to the conditional ones.
      tab = prop.table(tab, margin = seq(length(parents) + 1)[-1])

      structure(list(node = node, parents = parents, children = children,
        prob = tab), class = "bn.fit.dnode")

    }#FIT

  }#THEN
  else {

    fit = function(node) {

      # store the labels of the parents and the children to get
      # them only once.
      parents = x$nodes[[node]]$parents
      children = x$nodes[[node]]$children

      if (debug)
        cat("* fitting parameters of node", node, "(continuous).\n")

      if (length(parents) == 0) {

        # the only real parameters here are the intercept of the model
        # (which is the only coefficient) and the standard deviation of
        # the residuals.
        mean = mean(data[, node])
        coefs = c("(Intercept)" = mean)
        resid = data[, node] - mean
        sd = sd(data[, node])

        structure(list(node = node, parents = parents, children = children,
          coefficients = coefs, residuals = resid,
          fitted.values = rep(mean, n), sd = sd), class = "bn.fit.gnode")

      }#THEN
      else {

        if (debug)
          cat("  > found parents:", parents, "\n")

        # the relevant quantities of the normal distribution are the ones
        # usually found in an lm object: coefficients, fitted values and
        # residuals, plus the standard deviation of the latter.
        qr.x = qr(minimal.qr.matrix(data, parents))
        coefs = structure(qr.coef(qr.x, data[, node]), names = c("(Intercept)", parents))
        fitted = qr.fitted(qr.x, data[, node])
        resid = qr.resid(qr.x, data[, node])
        sd = sd(resid)

        structure(list(node = node, parents = parents, children = children,
          coefficients = c(coefs), residuals = resid,
          fitted.values = fitted, sd = sd), class = "bn.fit.gnode")

      }#ELSE

    }#FIT

  }#ELSE

  # fit the parameters of each node.
  fitted = sapply(names(x$nodes), fit, simplify = FALSE)
  # preserve any additional class of the original bn object.
  orig.class = class(x)
  class = c(orig.class[orig.class != "bn"], "bn.fit")
  # preserve the training node label from Bayesian network classifiers.
  if (x$learning$algo %in% classifiers)
    fitted = structure(fitted, class = class, training = x$learning$args$training)
  else
    fitted = structure(fitted, class = class)

  return(fitted)

}#BN.FIT.BACKEND

