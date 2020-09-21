
discretize.backend = function(data, method, breaks, ordered = FALSE, extra.args,
    debug = FALSE) {

  if (method %in% c("quantile", "interval")) {

    marginal.discretize.backend(data = data, method = method, breaks = breaks,
      ordered = ordered)

  }#THEN
  else if (method == "hartemink") {

    hartemink.discretization(data = data, breaks = breaks, ordered = ordered,
      initial.breaks = extra.args$ibreaks,
      initial.discretization = extra.args$idisc, debug = debug)

  }#ELSE

}#DISCRETIZE.BACKEND

hartemink.discretization = function(data, breaks, ordered, initial.breaks,
    initial.discretization, debug = FALSE) {

  # cache some useful quantities.
  nodes = names(data)
  nnodes = length(nodes)
  ndata = nrow(data)
  # check which type of data we are dealing with.
  type = data.type(data)
  # perform an initial discretization if needed.
  if (type == "continuous") {

    # make sure ordered is expanded, as expected by the backend implementing
    # marginal discretization methods.
    if (length(ordered) == 1)
      ordered = rep(ordered, ncol(data))

    discretized = discretize.backend(data = data, method = initial.discretization,
                    breaks = rep(initial.breaks, ncol(data)), extra.args = list(),
                    ordered = ordered, debug = FALSE)

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
      cur.levels = levels(.data.frame.column(discretized, node))

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
                    y = .data.frame.column(discretized, node),
                    ndata = ndata)[1]

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

dedup.backend = function(data, threshold, complete, debug = FALSE) {

  .Call(call_dedup,
        data = data,
        threshold = threshold,
        complete = complete,
        debug = debug)

}#DEDUP.BACKEND

marginal.discretize.backend = function(data, method, breaks, ordered = FALSE) {

  .Call(call_marginal_discretize,
        data = data,
        method = method,
        breaks = as.integer(breaks),
        ordered = ordered)

}#MARGINAL.DISCRETIZE.BACKEND

joint.discretize.backend = function(data, method, breaks, ordered = FALSE,
    initial.discretization = "quantile", initial.breaks) {

  .Call(call_joint_discretize,
        data = data,
        method = method,
        breaks = as.integer(breaks),
        ordered = ordered,
        initial.discretization = initial.discretization,
        initial.breaks = as.integer(initial.breaks))

}#JOINT.DISCRETIZE.BACKEND
