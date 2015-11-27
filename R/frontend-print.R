
# print method for class bn.
print.bn = function(x, ...) {

  params = names(x$learning$args)
  directed.arcs = length(which(which.directed(x$arcs)))
  undirected.arcs = (nrow(x$arcs) - directed.arcs)/2
  arcs = undirected.arcs + directed.arcs
  nodes = names(x$nodes)
  avg.mb = mean(sapply(nodes, function(n) { length(x$nodes[[n]]$mb) }))
  avg.nbr = mean(sapply(nodes, function(n) { length(x$nodes[[n]]$nbr) }))
  avg.ch = mean(sapply(nodes, function(n) { length(x$nodes[[n]]$children) }))

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (x$learning$algo %in% graph.generation.algorithms)
    cat("\n  Random/Generated Bayesian network\n\n")
  else if (x$learning$algo %in% constraint.based.algorithms)
    cat("\n  Bayesian network learned via Constraint-based methods\n\n")
  else if (x$learning$algo %in% score.based.algorithms)
    cat("\n  Bayesian network learned via Score-based methods\n\n")
  else if (x$learning$algo %in% hybrid.algorithms)
    cat("\n  Bayesian network learned via Hybrid methods\n\n")
  else if (x$learning$algo %in% mim.based.algorithms)
    cat("\n  Bayesian network learned via Pairwise Mutual Information methods\n\n")
  else if (x$learning$algo %in% classifiers)
    cat("\n  Bayesian network Classifier\n\n")
  else
    cat("\n  Bayesian network learned via [unknown] methods\n\n")

  # print the model string if possible, a short description otherwise.
  cat("  model:\n")

  if (undirected.arcs == 0)
    fcat(formula.backend(x))
  else if (directed.arcs == 0)
    cat("    [undirected graph]\n")
  else
    cat("    [partially directed graph]\n")

  wcat("  nodes:                                ", length(x$nodes))
  wcat("  arcs:                                 ", arcs)
  wcat("    undirected arcs:                    ", undirected.arcs)
  wcat("    directed arcs:                      ", directed.arcs)
  wcat("  average markov blanket size:          ", format(avg.mb, digits = 2, nsmall = 2))
  wcat("  average neighbourhood size:           ", format(avg.nbr, digits = 2, nsmall = 2))
  wcat("  average branching factor:             ", format(avg.ch, digits = 2, nsmall = 2))

  cat("\n")

  if (x$learning$algo %!in% graph.generation.algorithms) {

    wcat("  learning algorithm:                   ", method.labels[x$learning$algo])

    if (x$learning$algo %in% constraint.based.algorithms)
      wcat("  conditional independence test:        ", test.labels[x$learning$test])
    else if (x$learning$algo %in% score.based.algorithms)
      wcat("  score:                                ", score.labels[x$learning$test])
    else if (x$learning$algo %in% hybrid.algorithms) {

      if (x$learning$restrict %in% constraint.based.algorithms) {

        wcat("  constraint-based method:              ", method.labels[x$learning$restrict])
        wcat("  conditional independence test:        ", test.labels[x$learning$rstest])

      }#THEN
      else if (x$learning$restrict %in% mim.based.algorithms) {

        wcat("  pairwise mutual information method:   ", method.labels[x$learning$restrict])
        wcat("  mutual information estimator:         ", format(mi.estimator.labels[x$learning$args$estimator]))

      }#THEN

      wcat("  score-based method:                   ", method.labels[x$learning$maximize])
      wcat("  score:                                ", score.labels[x$learning$maxscore])

    }#THEN

    if ("prior" %in% params)
      wcat("  graph prior:                          ", prior.labels[x$learning$args$prior])
    if ("beta" %in% params)
      if (x$learning$args$prior != "cs")
        wcat("  beta sparsity parameter:              ", format(x$learning$args$beta))
      else
        wcat("  beta sparsity parameter:              ", "Completed Prior over Arcs")
    if ("alpha" %in% params)
      wcat("  alpha threshold:                      ", format(x$learning$args$alpha))
    if ("B" %in% params)
      wcat("  permutations:                         ", format(x$learning$args$B))
    if ("iss" %in% params)
      wcat("  imaginary sample size:                ", format(x$learning$args$iss))
    if ("phi" %in% params)
      wcat("  phi matrix structure:                 ", x$learning$args$phi)
    if ("k" %in% params)
      wcat("  penalization coefficient:             ", format(x$learning$args$k))
    if ("maxp" %in% params)
      wcat("  maximum parents:                      ", format(x$learning$args$maxp))

    if (x$learning$algo %in% c(mim.based.algorithms, classifiers)) {

      if ("estimator" %in% params)
        wcat("  mutual information estimator:         ", format(mi.estimator.labels[x$learning$args$estimator]))

      wcat("  training node:                        ", x$learning$args$training)

    }#THEN

    wcat("  tests used in the learning procedure: ", x$learning$ntests)

    if (!is.null(x$learning$optimized))
      wcat("  optimized:                            ", x$learning$optimized)

  }#THEN
  else {

    wcat("  generation algorithm:                 ", graph.generation.labels[x$learning$algo])

    if ("prob" %in% params)
      wcat("  arc sampling probability:             ", format(x$learning$args$prob))
    if ("burn.in" %in% params)
      wcat("  burn in length:                       ", format(x$learning$args$burn.in))
    if ("max.in.degree" %in% params)
      wcat("  maximum in-degree:                    ", format(x$learning$args$max.in.degree))
    if ("max.out.degree" %in% params)
      wcat("  maximum out-degree:                   ", format(x$learning$args$max.out.degree))
    if ("max.degree" %in% params)
      wcat("  maximum degree:                       ", format(x$learning$args$max.degree))
    if ("threshold" %in% params)
      wcat("  significance threshold:               ", format(x$learning$args$threshold))

  }#ELSE

  cat("\n")

  invisible(x)

}#PRINT.BN

# print method for class bn.fit.
print.bn.fit = function(x, order, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))
  # set the order in which the modes are to be printed.
  if (missing(order)) {

    # print the nodes in the order in which they are stored.
    order = seq_along(x)

  }#THEN
  else {

    # check the node ordering.
    check.nodes(order, graph = x, min.nodes = length(x), max.nodes = length(x))

  }#ELSE

  cat("\n  Bayesian network parameters\n")

  for (i in order)
    print(x[[i]])

  cat("\n")

  invisible(x)

}#PRINT.BN.FIT

# print method for class bn.fit.dnode.
print.bn.fit.dnode = function(x, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (is(x, "bn.fit.dnode"))
    cat("\n  Parameters of node", x$node, "(multinomial distribution)\n")
  else if (is(x, "bn.fit.onode"))
    cat("\n  Parameters of node", x$node, "(ordinal distribution)\n")

  cat("\nConditional probability table:\n", ifelse(length(x$parents) > 0, "\n", ""))
  print(x$prob)

  invisible(x)

}#PRINT.BN.FIT.DNODE

# print method for class bn.fit.onode.
print.bn.fit.onode = print.bn.fit.dnode

# print method for class bn.fit.gnode.
print.bn.fit.gnode = function(x, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  cat("\n  Parameters of node", x$node, "(Gaussian distribution)\n")

  cat("\nConditional density: ")
  if (length(x$parents) > 0)
    cat(paste(x$node, "|", paste(x$parents, sep = "", collapse = " + ")))
  else
    cat(x$node)

  cat("\nCoefficients:\n")
  print.default(format(c(x$coefficients)), print.gap = 2, right = TRUE, quote = FALSE)

  cat("Standard deviation of the residuals:", x$sd, "\n")

  invisible(x)

}#PRINT.BN.FIT.GNODE

# print method for class bn.fit.cgnode.
print.bn.fit.cgnode = function(x, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  cat("\n  Parameters of node", x$node, "(conditional Gaussian distribution)\n")

  cat("\nConditional density: ")
  if (length(x$parents) > 0)
    cat(paste(x$node, "|", paste(x$parents, sep = "", collapse = " + ")))
  else
    cat(x$node)

  cat("\nCoefficients:\n")
  print.default(format(x$coefficients), print.gap = 2, right = TRUE, quote = FALSE)

  cat("Standard deviation of the residuals:\n")
  print.default(format(x$sd), print.gap = 2, right = TRUE, quote = FALSE)

  cat("Discrete parents' configurations:\n")
  config = expand.grid(x$dlevels)
  rownames(config) = seq(from = 0, to = nrow(config) - 1)
  print.data.frame(config, print.gap = 2, right = TRUE, quote = FALSE)

  invisible(x)

}#PRINT.BN.FIT.CGNODE

# print method for a single cross-validation run.
print.bn.kcv = function(x, print.loss = TRUE, ...) {

  a = attributes(x)

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  cat("\n ", a$method,  "cross-validation for Bayesian networks\n\n")

  if (is.character(a$bn)) {

    wcat("  target learning algorithm:            ", method.labels[a$bn])

  }#THEN
  else {

    cat("  target network structure:\n")
    if (is(a$bn, "bn.naive"))
      cat("   [Naive Bayes Classifier]\n")
    else
      fcat(formula.backend(a$bn))

  }#ELSE

  wcat("  number of subsets:                    ", length(x))
  wcat("  loss function:                        ", loss.labels[a$loss])

  if ("target" %in% loss.extra.args[[a$loss]])
    wcat("  training node:                        ", a$args$target)

  if (print.loss) {

    wcat("  expected loss:                        ", format(a$mean))
    cat("\n")

  }#THEN

  invisible(x)

}#PRINT.BN.KCV

# print method for a list containing multiple cross-validation runs.
print.bn.kcv.list = function(x, ...) {

  losses = sapply(x, function(x) attr(x, "mean"))

  print.bn.kcv(x[[1]], print.loss = FALSE)
  wcat("  number of runs:                       ", length(x))
  wcat("  average loss over the runs:           ", format(mean(losses)))
  wcat("  standard deviation of the loss:       ", format(cgsd(losses)))
  cat("\n")

  invisible(x)

}#PRINT.BN.KCV.LIST
