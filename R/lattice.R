# lattice backend for plots aimed at gaussian bayesian networks.
lattice.discrete.backend = function(fitted, type, xlab, ylab, main, ...) {

  # check whether lattice is loaded, and try to load if it is not.
  if (!requireNamespace("lattice"))
    stop("this function requires the lattice package.")

  if (is(fitted, "bn.fit"))
    stop("only plots of single, discrete nodes are implemented.")
  if (!is(fitted, c("bn.fit.dnode", "bn.fit.onode")))
    stop("fitted must be an object of class 'bn.fit.dnode' or 'bn.fit.onode'.")

  if (type == "bar") {

    if (length(fitted$parents) == 0) {

      p = lattice::barchart(fitted$prob, xlab = xlab, ylab = ylab, main = main,
            panel = function(x, y, ...) {
              lattice::panel.grid(h = 0, v = -1)
              lattice::panel.barchart(x, y, ...)
            })

    }#THEN
    else {

      p = lattice::barchart(fitted$prob, groups = FALSE, as.table = TRUE,
            xlab = xlab, ylab = ylab, main = main,
            panel = function(x, y, ...) {
              lattice::panel.grid(h = 0, v = -1)
              lattice::panel.barchart(x, y, ...)
            })

    }#ELSE

  }#THEN
  else if (type == "dot") {

    if (length(fitted$parents) == 0) {

      p = lattice::dotplot(fitted$prob, xlab = xlab, ylab = ylab, main = main,
            type = c("p", "h"),
            panel = function(x, y, ...) {
              lattice::panel.grid(h = 0, v = -1)
              lattice::panel.dotplot(x, y, ...)
            })

    }#THEN
    else {

      p = lattice::dotplot(fitted$prob, groups = FALSE, as.table = TRUE,
            type = c("p", "h"), xlab = xlab, ylab = ylab, main = main,
            panel = function(x, y, ...) {
              lattice::panel.grid(h = 0, v = -1)
              lattice::panel.dotplot(x, y, ...)
            })

    }#ELSE

  }#THEN

  # print the plot explicitly, do not rely on auto-printing.
  print(p)

}#LATTICE.DISCRETE.BACKEND

# lattice backend for plots aimed at gaussian bayesian networks.
lattice.gaussian.backend = function(fitted, type, xlab, ylab, main, ...) {

  # check whether lattice is loaded, and try to load if it is not.
  if (!requireNamespace("lattice"))
    stop("this function requires lattice.")

  if (is(fitted, "bn.fit")) {

    # plot a panel for each node in the bayesian network.
    if (!is(fitted, "bn.fit.gnet"))
      stop("this plot is limited to Gaussian bayesian networks.")
    # check whether the residuals are there.
    if (any(!sapply(fitted, function(x) "residuals" %in% names(x))))
      stop("no residuals present in the bn.fit object.")

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

      p = lattice::qqmath(~ resid | node, data = temp,
            xlab = xlab, ylab = ylab, main = main,
            panel = function(x, ...) {
              lattice::panel.qqmathline(x, ...)
              lattice::panel.qqmath(x, ...)
            })

    }#THEN
    else if (type == "hist-dens") {

      p = lattice::histogram(~ resid | node, data = temp,
            xlab = xlab, ylab = ylab, main = main, type = "density",
            panel = function(x, ...) {
              lattice::panel.histogram(x, ...)
              lattice::panel.mathdensity(dmath = dnorm, col = "black",
                         args = list(mean = mean(x),
                                  sd = fitted[[lattice::panel.number()]]$sd))
            })

    }#THEN
    else if (type == "hist") {

      p = lattice::histogram( ~ resid | node, data = temp,
            xlab = xlab, ylab = ylab, main = main)

    }#THEN
    else if (type == "fitted") {

      # check whether the residuals are there.
      if (any(!sapply(fitted, function(x) "fitted.values" %in% names(x))))
        stop("no fitted values present in the bn.fit object.")

      p = lattice::xyplot(resid ~ fitted | node, data = temp,
            xlab = xlab, ylab = ylab, main = main,
            panel = function(x, ...) {
              lattice::panel.xyplot(x, ...)
              lattice::panel.abline(h = 0)
            })

    }#THEN

  }#THEN
  else if (is(fitted, c("bn.fit.gnode", "bn.fit.cgnode"))) {

    # check whether the residuals are there.
    if ("residuals" %!in% names(fitted))
      stop("no residuals present in the bn.fit.gnode object.")

    # print the equivalent plot for a single node.
    if (type == "qqplot") {

      f = formula(ifelse(is.null(fitted$configs),
            "~ residuals", "~ residuals | configs"))
      p = lattice::qqmath(f, data = fitted,
            xlab = xlab, ylab = ylab, main = main,
            panel = function(x, ...) {
              lattice::panel.qqmathline(x, ...)
              lattice::panel.qqmath(x, ...)
            })

    }#THEN
    else if (type == "hist-dens") {

      f = formula(ifelse(is.null(fitted$configs),
            "~ residuals", "~ residuals | configs"))
      p = lattice::histogram(f, data = fitted,
            xlab = xlab, ylab = ylab, main = main, type = "density",
            panel = function(x, ...) {
              lattice::panel.histogram(x, ...)
              lattice::panel.mathdensity(dmath = dnorm, col = "black",
                         args = list(mean = mean(x), sd = sd(x)))
            })

    }#THEN
    else if (type == "hist") {

      f = formula(ifelse(is.null(fitted$configs),
            "~ residuals", "~ residuals | configs"))
      p = lattice::histogram(f, data = fitted,
            xlab = xlab, ylab = ylab, main = main)

    }#THEN
    else if (type == "fitted") {

      # check whether the fitted values are there.
      if ("fitted.values" %!in% names(fitted))
        stop("no fitted values present in the bn.fit.gnode object.")

      f = formula(ifelse(is.null(fitted$configs),
            "residuals ~ fitted.values", "residuals ~ fitted.values | configs"))
      p = lattice::xyplot(f, data = fitted,
            xlab = xlab, ylab = ylab, main = main,
            panel = function(x, ...) {
              lattice::panel.xyplot(x, ...)
              lattice::panel.abline(h = 0)
            })

    }#THEN

  }#THEN
  else {

    stop("fitted must be an object of class 'bn.fit', 'bn.fit.gnode' and 'bn.fit.cgnode'.")

  }#ELSE

  # print the plot explicitly, do not rely on auto-printing.
  print(p)

}#LATTICE.GAUSSIAN.BACKEND

# lattice backend for cross-validation boxplots.
lattice.cv.bwplot = function(means, labels, losses, main, xlab, ylab,
    connect = FALSE) {

  # check whether lattice is loaded, and try to load if it is not.
  if (!requireNamespace("lattice"))
    stop("this function requires the lattice package.")
  
  # check the labels of the cross-validation objects, if any.
  if (!missing(xlab)) {

    if (!is.string.vector(xlab))
      stop("'xlab' must be a vector of character strings,",
           "the labels of the cross-validation objects.")
    if (length(xlab) != length(unique(labels)))
      stop("wrong number of labels for the cross-validation objects.")

  }#THEN
  else {

    xlab = as.character(seq(length(unique(labels))))

  }#ELSE

  # if no axis label is specified, use the loss function label when possible.
  if (missing(ylab)) {

    if (length(unique(losses)) == 1)
      ylab = loss.labels[losses[1]]
    else
      ylab = "loss"

  }#THEN

  # if no title is specified, leave it empty.
  if (missing(main))
    main = ""

  # convert the labels into a factor with the xlab as levels.
  labels = factor(labels, levels = unique(labels))
  levels(labels) = xlab

   p = lattice::bwplot(means ~ factor(labels),
         horizontal = FALSE, pch = 19, ylab = ylab, main = main,
         panel = function(...) {

           lattice::panel.grid(..., h = -1, v = 0, lty = 2)
           lattice::panel.bwplot(...)

           # connect the median points of the boxplots.
           if (connect)
             lattice::panel.average(..., fun = median, col.line = "black") 

         },
         par.settings = list(box.umbrella = list(lty = 1, col = "black"),
           box.rectangle = list(col = "black", fill = "ivory"),
           plot.symbol = list(col = "black")))

  # print the plot explicitly, do not rely on auto-printing.
  print(p)

}#LATTICE.CV.BWPLOT

