# lattice backend for plots aimed at gaussian bayesian networks.
lattice.discrete.backend = function(fitted, type, xlab, ylab, main, ...) {

  # check whether lattice is loaded, and try to load if it is not.
  if (!require(lattice))
    stop("this function requires the lattice package.")

  if (is(fitted, "bn.fit")) {

    # plot a panel for each node in the bayesian network.

    if (!is.fitted.discrete(fitted))
      stop("this plot is limited to discrete bayesian networks.")

    stop("only plots of single nodes are implemented.")

  }#THEN
  else if (class(fitted) == "bn.fit.dnode") {

    # print the equivalent plot for a single node.

    if (type == "bar") {

      if (length(fitted$parents) == 0) {

        barchart(fitted$prob, xlab = xlab, ylab = ylab, main = main,
          panel = function(x, y, ...) {
            panel.grid(h = 0, v = -1)
            panel.barchart(x, y, ...)
          })

      }#THEN
      else {

        barchart(fitted$prob, groups = FALSE, as.table = TRUE,
          xlab = xlab, ylab = ylab, main = main,
          panel = function(x, y, ...) {
            panel.grid(h = 0, v = -1)
            panel.barchart(x, y, ...)
          })

      }#ELSE

    }#THEN
    else if (type == "dot") {

      if (length(fitted$parents) == 0) {

        dotplot(fitted$prob, xlab = xlab, ylab = ylab, main = main,
          type = c("p", "h"),
          panel = function(x, y, ...) {
            panel.grid(h = 0, v = -1)
            panel.dotplot(x, y, ...)
          })

      }#THEN
      else {

        dotplot(fitted$prob, groups = FALSE, as.table = TRUE,
          type = c("p", "h"),
          xlab = xlab, ylab = ylab, main = main,
          panel = function(x, y, ...) {
            panel.grid(h = 0, v = -1)
            panel.dotplot(x, y, ...)
          })

      }#ELSE

    }#THEN

  }#THEN
  else {

    stop("fitted must be an object of class 'bn.fit.dnode'.")

  }#ELSE

}#LATTICE.DISCRETE.BACKEND

# lattice backend for plots aimed at gaussian bayesian networks.
lattice.gaussian.backend = function(fitted, type, xlab, ylab, main, ...) {

  # check whether lattice is loaded, and try to load if it is not.
  if (!require(lattice))
    stop("this function requires lattice.")

  if (is(fitted, "bn.fit")) {

    # plot a panel for each node in the bayesian network.

    if (!is.fitted.continuous(fitted))
      stop("this plot is limited to Gaussian bayesian networks.")

    nodes = names(fitted)
    nrows = length(fitted[[1]]$residuals)

    if (type %in% c("fitted", "fitted-std")) {

      temp = data.frame(
        resid = unlist(lapply(fitted, "[[", "residuals" )),
        fitted = unlist(lapply(fitted, "[[", "fitted.values" )),
        node = unlist(lapply(nodes, rep, times = nrows))
      )

    }#THEN
    else {

      temp = data.frame(
        resid = unlist(lapply(fitted, "[[", "residuals" )),
        node = unlist(lapply(nodes, rep, times = nrows))
      )

    }#ELSE

    if (type == "qqplot") {

      qqmath(~ resid | node, data = temp,
        xlab = xlab, ylab = ylab, main = main,
        panel = function(x, ...) {
          panel.qqmathline(x, ...)
          panel.qqmath(x, ...)
        })

    }#THEN
    else if (type == "hist-dens") {

      histogram(~ resid | node, data = temp,
        xlab = xlab, ylab = ylab, main = main, type = "density",
        panel = function(x, ...) {
          panel.histogram(x, ...)
          panel.mathdensity(dmath = dnorm, col = "black",
            args = list(mean = mean(x),
              sd = fitted[[panel.number()]]$sd))
        })

    }#THEN
    else if (type == "hist") {

      histogram( ~ resid | node, data = temp,
        xlab = xlab, ylab = ylab, main = main)

    }#THEN
    else if (type == "fitted") {

      xyplot(resid ~ fitted | node, data = temp,
        xlab = xlab, ylab = ylab, main = main,
        panel = function(x, ...) {
          panel.xyplot(x, ...)
          panel.abline(h = 0)
        })

    }#THEN

  }#THEN
  else if (class(fitted) == "bn.fit.gnode") {

    # print the equivalent plot for a single node.

    if (type == "qqplot") {

      qqmath(~ residuals , data = fitted,
        xlab = xlab, ylab = ylab, main = main,
        panel = function(x, ...) {
          panel.qqmathline(x, ...)
          panel.qqmath(x, ...)
        })

    }#THEN
    else if (type == "hist-dens") {

      histogram(~ residuals, data = fitted,
        xlab = xlab, ylab = ylab, main = main, type = "density",
        panel = function(x, ...) {
          panel.histogram(x, ...)
          panel.mathdensity(dmath = dnorm, col = "black",
            args = list(mean = mean(x), sd = sd(x)))
        })

    }#THEN
    else if (type == "hist") {

      histogram(~ residuals, data = fitted,
        xlab = xlab, ylab = ylab, main = main)

    }#THEN
    else if (type == "fitted") {

      xyplot(residuals ~ fitted.values, data = fitted,
        xlab = xlab, ylab = ylab, main = main,
        panel = function(x, ...) {
          panel.xyplot(x, ...)
          panel.abline(h = 0)
        })

    }#THEN

  }#THEN
  else {

    stop("fitted must be an object of class 'bn.fit' or 'bn.fit.gnode'.")

  }#ELSE

}#LATTICE.GAUSSIAN.BACKEND

