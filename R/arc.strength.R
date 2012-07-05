
# compute arcs' strength as the p-value of the test for their removal.
arc.strength.test = function(network, data, test, alpha, B, debug = FALSE) {

  drop = function(arc) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* computing strength for arc", arc[1], "->", arc[2], ".\n")

    }#THEN

    parents =
      network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    a = conditional.test(arc[1], arc[2], parents, data = data, test = test,
          B = B, alpha = alpha)

    if (debug) {

      cat("  > testing", arc[1], "->", arc[2],
        "with conditioning set '", parents, "'.\n")
      cat("    > p-value is", a, ".\n")

    }#THEN

    return(a)

  }#DROP

  if (debug) {

    cat("----------------------------------------------------------------\n")
    print(network)

  }#THEN

  # populate the strength data frame.
  strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
               stringsAsFactors = FALSE)

  return(strength)

}#ARC.STRENGTH.TEST

# compute arcs' strength as the difference in score caused by their removal.
arc.strength.score = function(network, data, score, extra, debug = FALSE) {

  drop = function(arc) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* computing strength for arc", arc[1], "->", arc[2], ".\n")

    }#THEN

    better = score.delta(arc = arc, network = network, data = data,
               score = score, score.delta = 0,
               reference.score = reference.score, op = "drop",
               extra = extra)

    if (debug) {

      cat("  > updated score for node", arc[2], "is", better$updates, ".\n")
      cat("  > score delta", better$delta, ".\n")

    }#THEN

    return(better$delta)

  }#DROP

  # cache nodes' labels.
  nodes = names(data)
  # set the reference score.
  reference.score = per.node.score(network = network, score = score,
                      nodes = nodes, extra.args = extra, data = data)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    print(network)
    cat("* current score:", sum(reference.score), "\n")
    cat("----------------------------------------------------------------\n")
    cat("* original scores of the nodes of the graphs are:\n")
    for (n in nodes)
      cat("  > original score for node", n, "is", reference.score[n], ".\n")

  }#THEN

  # populate the strength data frame.
  strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
               stringsAsFactors = FALSE)

  return(strength)

}#ARC.STRENGTH.SCORE

# compute arcs' strength as the relative frequency in boostrapped networks.
arc.strength.boot = function(data, R, m, algorithm, algorithm.args, arcs,
    cpdag, debug = FALSE) {

  # allocate and initialize an empty adjacency matrix.
  prob = matrix(0, ncol = ncol(data), nrow = ncol(data))
  # get the names of the variables in the data set.
  nodes = names(data)

  for (r in seq_len(R)) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* bootstrap replicate", r, ".\n")

    }#THEN

    # perform the resampling with replacement.
    resampling = sample(nrow(data), m, replace = TRUE)

    # user-provided lists of manipulated observations for the mbde score must
    # be remapped to match the bootstrap sample.
    if ((algorithm.args$score == "mbde") && !is.null(algorithm.args$exp)) {

      algorithm.args$exp = lapply(algorithm.args$exp, function(x) {

        x = match(x, resampling)
        x = x[!is.na(x)]

      })

    }#THEN

    # generate the r-th bootstrap sample.
    replicate = data[resampling, , drop = FALSE]

    # learn the network structure from the bootstrap sample.
    net = do.call(algorithm, c(list(x = replicate), algorithm.args))

    # switch to the equivalence class if asked to.
    if (cpdag)
      net = cpdag(net)

    if (debug) {

      cat("* learning bayesian network structure.\n")
      print(net)

    }#THEN

    # update the counters in the matrix: undirected arcs are counted half
    # for each direction, so that when summing up strength and direction
    # they get counted once instead of twice.
    # BEWARE: in-place modification of prob!
    .Call("bootstrap_strength_counters",
          prob = prob,
          weight = 1,
          arcs = net$arcs,
          nodes = nodes,
          PACKAGE = "bnlearn")

  }#FOR

  # rescale the counters to probabilities.
  prob = prob / R

  res = .Call("bootstrap_arc_coefficients",
              prob = prob,
              nodes = nodes,
              PACKAGE = "bnlearn")

  if (!is.null(arcs))
    return(match.arcs.and.strengths(arcs, nodes, res, keep = TRUE))
  else
    return(res)

}#ARC.STRENGTH.BOOT

# compute arcs' strength as the relative frequency in a custom list of
# network structures/arc sets.
arc.strength.custom = function(custom.list, nodes, arcs, cpdag, weights = NULL,
    debug = FALSE) {

  # allocate and initialize an empty adjacency matrix.
  prob = matrix(0, ncol = length(nodes), nrow = length(nodes))
  # get the number of networks to average.
  R = length(custom.list)
  # set the total weight mass.
  weight.sum = sum(weights)

  for (r in seq_len(R)) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* considering element", r, ".\n")

    }#THEN

    # get the arc set on way or the other.
    if (is(custom.list[[r]], "bn"))
      net.arcs = custom.list[[r]]$arcs
    else
      net.arcs = custom.list[[r]]

    # switch to the equivalence class if asked to.
    if (cpdag)
      net.arcs = cpdag.arc.backend(nodes, net.arcs)

    # update the counters in the matrix: undirected arcs are counted half
    # for each direction, so that when summing up strength and direction
    # they get counted once instead of twice.
    # BEWARE: in-place modification of prob!
    .Call("bootstrap_strength_counters",
          prob = prob,
          weight = weights[r],
          arcs = net.arcs,
          nodes = nodes,
          PACKAGE = "bnlearn")

  }#FOR

  # rescale the counters to probabilities.
  prob = prob / weight.sum

  res = .Call("bootstrap_arc_coefficients",
              prob = prob,
              nodes = nodes,
              PACKAGE = "bnlearn")

  if (!is.null(arcs))
    return(match.arcs.and.strengths(arcs, nodes, res, keep = TRUE))
  else
    return(res)

}#ARC.STRENGTH.CUSTOM

# match arc sets and strength coefficients.
match.arcs.and.strengths = function(arcs, nodes, strengths, keep = FALSE) {

  if (nrow(strengths) < nrow(arcs))
    stop("insufficient number of strength coefficients.")

  a_hash = interaction(arcs[, "from"], arcs[, "to"])
  s_hash = interaction(strengths[, "from"], strengths[, "to"])

  if (keep) {

    s = strengths[match(a_hash, s_hash), , drop = FALSE]
    coef = s$strength

  }#THEN
  else {

    s = strengths[match(a_hash, s_hash), "strength"]
    coef = s

  }#ELSE

  if (any(is.na(coef))) {

    missing = apply(arcs[is.na(coef), , drop = FALSE], 1,
                function(x) { paste(" (", x[1], ", ", x[2], ")", sep = "")  })

    stop("the following arcs do not have a corresponding strength coefficients:",
      missing, ".")

  }#THEN

  return(s)

}#MATCH.ARCS.AND.STRENGTH

# convert an arc strength object to the corresponding line widths for plotting.
strength2lwd = function(strength, threshold, cutpoints, mode, arcs = NULL,
    debug = FALSE) {

  # sanitize user-defined cut points, if any.
  if (!missing(cutpoints)) {

    if (!is.numeric(cutpoints) || any(is.nan(cutpoints)))
      stop("cut points must be numerical values.")
    if (length(strength) <= length(cutpoints))
      stop("there are at least as many cut points as strength values.")

  }#THEN

  if (debug) {

    cat("* using threshold:", threshold, "\n")
    cat("* reported arc strength are:\n")
    if (!is.null(arcs))
      print(cbind(arcs, strength))
    else
      print(strength)

  }#THEN

  if (mode %in% c("test", "bootstrap")) {

    # bootstrap probabilities work like p-values, only reversed.
    if (mode == "bootstrap") {

      strength = 1 - strength
      threshold = 1 - threshold

    }#THEN

    # use user-defined cut points if available.
    if (missing(cutpoints))
      cutpoints = unique(c(0, threshold/c(10, 5, 2, 1.5, 1), 1))
    else
      cutpoints = sort(cutpoints)

    # p-values are already scaled, so the raw quantiles are good cut points.
    arc.weights = cut(strength, cutpoints, labels = FALSE, include.lowest = TRUE)

    arc.weights = length(cutpoints) - arc.weights

  }#THEN
  else if (mode == "score") {

    # score deltas are defined on a reversed scale (the more negative
    # the better); change their sign (and that of the threshold)
    # for simplicity
    threshold = -threshold
    strength = -strength

    # define a set of cut points using the quantiles from the empirical
    # distribution of the negative score deltas (that is, the ones
    # corresponding to significant arcs) or use user-defined ones.
    if (missing(cutpoints)) {

      significant = strength[strength > threshold]
      q = quantile(significant, c(0.50, 0.75, 0.90, 0.95, 1), names = FALSE)
      cutpoints = sort(c(-Inf, threshold, unique(q), Inf))

    }#THEN

    arc.weights = cut(strength, cutpoints, labels = FALSE)

  }#THEN

  # arcs beyond the significance threshold are given a negative weight,
  # so that graphviz.backend() will draw them as dashed lines.
  arc.weights[arc.weights == 1] = -1

  if (debug) {

    cat("* using cut points for strength intervals:\n")
    if (mode == "boot")
      print(1 - cutpoints)
    else
      print(cutpoints)
    cat("* arc weights:", arc.weights, "\n")

  }#THEN

  return(arc.weights)

}#STRENGTH2LWD

threshold = function(strength, method = "l1") {

  e = ecdf(strength$strength)
  u = knots(e)

  if (method == "l1") {

    norm = function(p)
      sum( diff(unique(c(0, u, 1))) * abs(e(unique(c(0, u[u < 1]))) - p))

  }#THEN

  p0 = optimize(f = norm, interval = c(0, 1))$minimum

  # double check the boundaries, they are legal solutions but optimize() does
  # not check them.
  if (norm(1) < norm(p0))
    p0 = 1
  if (norm(0) < norm(p0))
    p0 = 0

  quantile(strength$strength, p0, type = 1, names = FALSE)

}#THRESHOLD

