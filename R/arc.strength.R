
# compute arcs' strength as the p-value of the test for their removal.
arc.strength.test = function(network, data, test, alpha, B, complete,
    debug = FALSE) {

  drop = function(arc) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* computing strength for arc", arc[1], "->", arc[2], ".\n")

    }#THEN

    parents =
      network$nodes[[arc[2]]]$parents[network$nodes[[arc[2]]]$parents != arc[1]]

    a = indep.test(arc[1], arc[2], parents, data = data, test = test,
          B = B, alpha = alpha, complete = complete)

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
  if (nrow(network$arcs) == 0) {

    strength = as.data.frame(list(from = character(0), to = character(0),
                 strength = numeric(0)), stringsAsFactors = FALSE)

  }#THEN
  else {

    strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
                 stringsAsFactors = FALSE)

  }#ELSE

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
               extra = extra, decomposable = decomp)

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
                      targets = nodes, extra.args = extra, data = data)
  # check whether the score is decomposable.
  decomp = is.score.decomposable(score, extra)

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
  if (nrow(network$arcs) == 0) {

    strength = as.data.frame(list(from = character(0), to = character(0),
                 strength = numeric(0)), stringsAsFactors = FALSE)

  }#THEN
  else {

    strength = data.frame(network$arcs, strength = apply(network$arcs, 1, drop),
                 stringsAsFactors = FALSE)

  }#ELSE

  return(strength)

}#ARC.STRENGTH.SCORE

# compute arcs' strength as the relative frequency in boostrapped networks.
arc.strength.boot = function(data, cluster = NULL, R, m, algorithm,
    algorithm.args, arcs, cpdag, debug = FALSE) {

  bootstrap.batch = function(data, R, m, algorithm, algorithm.args, arcs,
    cpdag, debug) {

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

        algorithm.args$exp =
          lapply(algorithm.args$exp,
            function(x) which(resampling %in% x))

      }#THEN

      # generate the r-th bootstrap sample.
      replicate = data[resampling, , drop = FALSE]

      # learn the network structure from the bootstrap sample.
      net = do.call(algorithm, c(list(x = replicate), algorithm.args))

      # switch to the equivalence class if asked to, but preserving whitelisted
      # and blacklisted arcs.
      if (cpdag)
        net = cpdag.backend(net, wlbl = TRUE, debug = FALSE)

      if (debug) {

        cat("* learning bayesian network structure.\n")
        print(net)

      }#THEN

      # update the counters in the matrix: undirected arcs are counted half
      # for each direction, so that when summing up strength and direction
      # they get counted once instead of twice.
      # BEWARE: in-place modification of prob!
      .Call(call_bootstrap_strength_counters,
            prob = prob,
            weight = 1,
            arcs = net$arcs,
            nodes = nodes)

    }#FOR

    # rescale the counters to probabilities.
    prob = prob / R

    res = .Call(call_bootstrap_arc_coefficients,
                prob = prob,
                nodes = nodes)

    if (!is.null(arcs))
      return(match.arcs.and.strengths(arcs, nodes, res, keep = TRUE))
    else
      return(res)

  }#BOOTSTRAP.BATCH

  if (!is.null(cluster)) {

    # get the number of slaves.
    s = nSlaves(cluster)

    res = parallel::parLapplyLB(cluster, rep(ceiling(R / s), s),
            bootstrap.batch, data = data, m = m, arcs = arcs,
            algorithm = algorithm, algorithm.args = algorithm.args,
            cpdag = cpdag, debug = debug)

    .Call(call_bootstrap_reduce,
          x = res)

  }#THEN
  else {

    bootstrap.batch(data = data, R = R, m = m, algorithm = algorithm,
      algorithm.args = algorithm.args, arcs = arcs, cpdag = cpdag,
      debug = debug)

  }#THEN

}#ARC.STRENGTH.BOOT

# compute an approximation of arc and direction strength from the Bayes factors
# that can be computed from a single MAP network.
bf.strength.backend = function(x, data, score, extra.args, precBits = 200,
  debug = FALSE) {

  # construct all pairs of nodes.
  nodes = names(x$nodes)
  all.pairs = structure(subsets(nodes, 2), dimnames = list(NULL, c("from", "to")))
  # prepare the return value.
  all.arcs = expand.grid(to = nodes, from = nodes, stringsAsFactors = FALSE)
  all.arcs = all.arcs[all.arcs$from != all.arcs$to, ]
  all.arcs = data.frame(all.arcs[, 2:1], strength = numeric(nrow(all.arcs)),
               direction = numeric(nrow(all.arcs)))

  # if a score is not decomposable, we stil assume it is possible to compute
  # score deltas by using only the two incident nodes as targets in
  # per.node.score().
  decomposable = is.score.decomposable(score, extra.args)
  amat = amat(x)

  if (debug) {

    cat("----------------------------------------------------------------\n")
    print(x)
    cat("----------------------------------------------------------------\n")

  }#THEN

  for (i in seq(nrow(all.pairs))) {

    arc = all.pairs[i, ]

    if (debug)
      cat("* considering the arc between", arc["from"], "and", arc["to"], ".\n")

    # first compute the score of the arc without the arc; it is always possible
    # to do that.
    dag.base = drop.arc(x, arc["from"], arc["to"])
    score.base = Rmpfr::mpfr(per.node.score(dag.base, data = data, score = score,
                 targets = arc, extra.args = extra.args), precBits = precBits)

    if (debug) {

      cat("  > no arc between", arc["from"], "and", arc["to"],
        ": base scores =", Rmpfr::asNumeric(score.base), ".\n")

    }#THEN

    try.arc = function(arc, score.base) {

      if (!has.path(arc[2], arc[1], nodes, amat, exclude.direct = TRUE)) {

        new.dag = set.arc(x, arc[1], arc[2], check.cycles = FALSE)

        if (decomposable) {

          new.score = Rmpfr::mpfr(per.node.score(new.dag, data = data,
                       score = score, targets = arc[2], extra.args = extra.args),
                       precBits = precBits)
          bf = exp(sum(new.score - score.base[2]))

        }#THEN
        else {

          new.score = Rmpfr::mpfr(per.node.score(new.dag, data = data,
                       score = score, targets = arc, extra.args = extra.args),
                       precBits = precBits)
          bf = exp(sum(new.score - score.base))

        }#ELSE

        if (debug) {

          cat("  > arc", arc[1], "->", arc[2], ": new score(s) =",
            Rmpfr::asNumeric(new.score), ", BF =",
            Rmpfr::format(bf, digits = 4), ".\n")

        }#THEN

      }#THEN
      else {

        bf = Rmpfr::mpfr(0, precBits = precBits)

        if (debug)
          cat("  > arc", arc[1], "->", arc[2], ": cycles!\n")

      }#ELSE

      return(bf)

    }#TRY.ARC

    # try to add the arc in one direction (arbitrarily called "forward").
    fw = try.arc(arc, score.base = score.base)
    # try to add the arc in the other direction (arbitrarily called "backward").
    bw = try.arc(rev(arc), score.base = rev(score.base))

    num = c(bw = bw, no = 1, fw = fw)
    norm = 1 + bw + fw

    # handle corner cases in which unnormalized probabilities are not finite;
    # probability is distributed uniformly over those that are infinite.
    if (!is.finite(norm)) {

      finite = is.finite(num)
      num[finite] = 0
      num[!finite] = 1/length(which(!finite))
      norm = 1

    }#THEN

    # normalize and save arc strength and directions.
    normalized = num / norm
    strength = Rmpfr::asNumeric(normalized[1] + normalized[3])
    if (normalized[2] == 1)
      direction = 0
    else
      direction = Rmpfr::asNumeric(normalized[3] / (normalized[1] + normalized[3]))

    all.arcs[(all.arcs$from == arc["from"]) & (all.arcs$to == arc["to"]),
      c("strength", "direction")] = c(strength, direction)
    all.arcs[(all.arcs$from == arc["to"]) & (all.arcs$to == arc["from"]),
      c("strength", "direction")] = c(strength, 1 - direction)

  }#FOR

  return(structure(all.arcs, row.names = seq(nrow(all.arcs))))

}#BF.STRENGTH.BACKEND

# compute arcs' strength as the relative frequency in a custom list of
# network structures/arc sets.
arc.strength.custom = function(custom.list, nodes, arcs, cpdag, weights = NULL,
    illegal = NULL, debug = FALSE) {

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

    # get the network structure one way or the other.
    if (is(custom.list[[r]], "bn"))
      net = custom.list[[r]]
    else if (is(custom.list[[r]], "bn.fit"))
      net = bn.net(custom.list[[r]])
    else {

      net = empty.graph(nodes)
      arcs(net) = custom.list[[r]]

    }#THEN

    # switch to the equivalence class if asked to.
    if (cpdag)
      net = cpdag.backend(net, moral = FALSE, wlbl = TRUE)

    # update the counters in the matrix: undirected arcs are counted half
    # for each direction, so that when summing up strength and direction
    # they get counted once instead of twice.
    # BEWARE: in-place modification of prob!
    .Call(call_bootstrap_strength_counters,
          prob = prob,
          weight = weights[r],
          arcs = net$arcs,
          nodes = nodes)

  }#FOR

  # rescale the counters to probabilities.
  prob = prob / weight.sum

  res = .Call(call_bootstrap_arc_coefficients,
              prob = prob,
              nodes = nodes)

  if (!is.null(arcs))
    return(match.arcs.and.strengths(arcs, nodes, res, keep = TRUE))
  else
    return(res)

}#ARC.STRENGTH.CUSTOM

# compute the significance threshold for Friedman's confidence.
threshold = function(strength, method = "l1") {

  # do not blow up with graphs with only 1 node.
  if (nrow(strength) == 0)
    return(0)

  e = ecdf(strength$strength)
  u = knots(e)

  if (method == "l1") {

    norm = function(p)
      sum( diff(unique(c(0, u, 1))) * abs(e(unique(c(0, u[u < 1]))) - p))

  }#THEN

  p0 = optimize(f = norm, interval = c(0, 1))$minimum

  # double-check the boundaries, they are legal solutions but optimize() does
  # not check them.
  if (norm(1) < norm(p0))
    p0 = 1
  if (norm(0) < norm(p0))
    p0 = 0

  quantile(strength$strength, p0, type = 1, names = FALSE)

}#THRESHOLD

# weighted average of bn.strength objects.
mean.strength = function(strength, nodes, weights) {

  .Call(call_mean_strength,
        strength = strength,
        nodes = nodes,
        weights = weights)

}#MEAN.STRENGTH
