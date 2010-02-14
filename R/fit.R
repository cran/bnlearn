
# fit the parameters of the bayesian network for a given network stucture.
bn.fit.backend = function(x, data, debug) {

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

      # the parameters of the multinomial distribution are the probabilities
      # of the levels of the node and the configurations of its parents.
      tab = table(data[, c(node, parents), drop = FALSE])
      # switch from the joint probabilities to the conditional ones.
      tab = prop.table((tab), margin = seq(length(parents) + 1)[-1])

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
        intercept = rep(1, nrow(data))
        qr.x = qr(cbind(intercept, data[, parents]))
        coefs = structure(qr.coef(qr.x, data[, node]), names = c("(Intercept)", parents))
        fitted = qr.fitted(qr.x, data[, node])
        resid = qr.resid(qr.x, data[, node])
        sd = sd(resid)

        structure(list(node = node, parents = parents, children = children,
          coefficients = coefs, residuals = resid, 
          fitted.values = fitted, sd = sd), class = "bn.fit.gnode")


      }#ELSE

    }#FIT

  }#ELSE

  # fit the parameters of each node.
  structure(sapply(names(x$nodes), fit, simplify = FALSE), class = "bn.fit")

}#BN.FIT.BACKEND
