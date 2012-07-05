
discretize.backend = function(data, method, breaks, extra.args, debug = FALSE) {

  if (method == "quantile")
    quantile.discretization(data = data, breaks = breaks)
  else if (method == "interval")
    interval.discretization(data = data, breaks = breaks)
  else if (method == "hartemink")
    hartemink.discretization(data = data, breaks = breaks,
      initial.breaks = extra.args$ibreaks,
      initial.discretization = extra.args$idisc, debug = debug)

}#DISCRETIZE.BACKEND

quantile.discretization = function(data, breaks) {

  discretized = lapply(seq(ncol(data)), function(x) {

    breaks = breaks[x]
    y = minimal.data.frame.column(data, x)

    # do not touch discrete variables.
    if (is(y, "factor"))
      return(y)
    # compute the quantiles for the variable.
    quantiles = quantile(y, probs = seq(from = 0, to = breaks)/breaks)
    # check whther the quantiles are unique.
    if (any(duplicated(quantiles)))
      stop("unable to discretize ", names(data)[x], " in ", breaks,
           " intervals, some quantiles are not unique.")
    # cut the range using the quantiles as break points.
    cut(y, breaks = quantiles, include.lowest = TRUE)

  })
  # convert the return value to a data frame.
  discretized = minimal.data.frame(discretized)
  names(discretized) = names(data)

  return(discretized)

}#QUANTILE.DISCRETIZATION

interval.discretization = function(data, breaks) {

  # cut the range into intervals and assign observations.
  discretized = lapply(seq(ncol(data)), function(x) {

    breaks = breaks[x]
    y = minimal.data.frame.column(data, x)

    # do not touch discrete variables.
    if (is(y, "factor"))
      return(y)
    # cut the range with the given number of break points.
    cut(y, breaks = breaks, include.lowest = TRUE)

  })
  # convert the return value to a data frame.
  discretized = minimal.data.frame(discretized)
  names(discretized) = names(data)

  return(discretized)

}#INTERVAL.DISCRETIZATION

hartemink.discretization = function(data, breaks, initial.breaks,
    initial.discretization, debug = FALSE) {

  # cache some useful quantities.
  nodes = names(data)
  nnodes = length(nodes)
  ndata = nrow(data)
  # perform an initial discretization if needed.
  if (is.data.continuous(data)) {

    discretized = discretize.backend(data = data, method = initial.discretization,
                    breaks = rep(initial.breaks, ncol(data)), extra.args = list(),
                    debug = FALSE)

  }#THEN
  else {

    discretized = data

  }#THEN

  # nothing to do, move along.
  if (initial.breaks == breaks)
    return(discretized)

  # conting down to the desired number of levels ...
  for (nlevels in seq(from = initial.breaks, to = breaks + 1)) {

    if (debug)
      cat("* considering", nlevels, "levels.\n")

    # ... for each variable ...
    for (node in nodes) {

      if (debug)
        cat("* considering variable", node, ".\n")

      total.mutual.information = numeric(nnodes - 1)
      cur.levels = levels(minimal.data.frame.column(discretized, node))

      # ... for each pair of levels ...
      for (collapsing in seq(nlevels - 1)) {


        if (debug)
          cat("  > collapsing levels", cur.levels[collapsing], "and", cur.levels[collapsing + 1], ".\n")

        # ... collapse them ...
        collapsed = discretized[, node]
        collapsed[collapsed == cur.levels[collapsing + 1]] = cur.levels[collapsing]

        # ... and compute the total mutual information.
	total.mutual.information[collapsing] = sum(sapply(nodes[nodes != node],
          function(node) {

            mi.test(x = collapsed,
                    y = minimal.data.frame.column(discretized, node),
                    ndata = ndata)

          }))

        if (debug)
          cat("    > total mutual information is", total.mutual.information[collapsing], ".\n")

      }#FOR

      # get which collapsing gives the best restults.
      bestop = which.max(total.mutual.information)

      if (debug)
        cat("@ best collapsing is", bestop, ".\n")

      # build the new set of levels.
      from = strsplit(cur.levels[bestop], "\\(|\\[|,|\\]|\\)")[[1]][2]
      to = strsplit(cur.levels[bestop + 1], "\\(|\\[|,|\\]|\\)")[[1]][3]
      new.levels = cur.levels

      if (is.na(from) || is.na(to)) {

        cur.levels[bestop] = gsub("\\[|\\]", "", cur.levels[bestop])
        cur.levels[bestop + 1] = gsub("\\[|\\]", "", cur.levels[bestop + 1])
        new.levels[bestop] = new.levels[bestop + 1] = paste("[", cur.levels[bestop], ":", cur.levels[bestop + 1], "]", sep = "")

      }#THEN
      else {

        new.levels[bestop] = new.levels[bestop + 1] = paste("(", from, ",", to, "]", sep = "")

      }#ELSE

      # collapse the levels in the discretized data.
      levels(discretized[, node]) = new.levels

    }#FOR

  }#FOR

  return(discretized)

}#HARTEMINK.DISCRETIZATION

