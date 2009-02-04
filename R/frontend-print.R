
# print method for class bn.
print.bn = function(x, ...) {

  params = names(x$learning$args)
  directed.arcs = length(which(!which.undirected(x$arcs)))
  undirected.arcs = (nrow(x$arcs) - directed.arcs)/2
  arcs = undirected.arcs + directed.arcs
  avg.mb = mean(sapply(nodes(x), function(n) { length(x$nodes[[n]]$mb) }))
  avg.nbr = mean(sapply(nodes(x), function(n) { length(x$nodes[[n]]$nbr) }))
  avg.ch = mean(sapply(nodes(x), function(n) { length(x$nodes[[n]]$children) }))

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  if (x$learning$test == "none")
    cat("\n  Randomly generated Bayesian network\n\n")
  else if (x$learning$test %in% names(test.labels))
    cat("\n  Bayesian network learned via Conditional Independence methods\n\n")
  else
    cat("\n  Bayesian network learned via Scoring methods\n\n")

  cat("  model:\n   ", ifelse(is.dag(x$arcs, names(x$nodes)), formula.backend(x),
      "[partially directed graph]"), "\n")

  cat("  nodes:                                ", length(x$nodes), "\n")
  cat("  arcs:                                 ", arcs, "\n")
  cat("    undirected arcs:                    ", undirected.arcs, "\n")
  cat("    directed arcs:                      ", directed.arcs, "\n")
  cat("  average markov blanket size:          ", format(avg.mb, digits = 2, nsmall = 2), "\n")
  cat("  average neighbourhood size:           ", format(avg.nbr, digits = 2, nsmall = 2), "\n")
  cat("  average branching factor:             ", format(avg.ch, digits = 2, nsmall = 2), "\n")

  cat("\n")

  if (x$learning$test != "none") {

    cat("  learning algorithm:                   ", method.labels[x$learning$algo], "\n")
  
    if (x$learning$test %in% names(test.labels))
      cat("  conditional independence test:        ", test.labels[x$learning$test], "\n")
    else
      cat("  score:                                ", score.labels[x$learning$test], "\n")

    if ("alpha" %in% params)
      cat("  alpha threshold:                      ", x$learning$args$alpha, "\n")
    if ("iss" %in% params)
      cat("  imaginary sample size:                ", x$learning$args$iss, "\n")
    if ("phi" %in% params)
      cat("  phi matrix structure:                 ", x$learning$args$phi, "\n")
    if ("k" %in% params)
      cat("  penalization coefficient:             ", x$learning$args$k, "\n")

    cat("  tests used in the learning procedure: ", x$learning$ntests, "\n")

  }#THEN
  else {

    cat("  generation algorithm:                 ", graph.generation.labels[x$learning$algo], "\n")

    if ("prob" %in% params)
      cat("  arc sampling probability:             ", x$learning$args$prob, "\n")
    if ("burn.in" %in% params)
      cat("  burn in length:                       ", x$learning$args$burn.in, "\n")
    if ("max.in.degree" %in% params)
      cat("  maximum in-degree:                    ", x$learning$args$max.in.degree, "\n")
    if ("max.out.degree" %in% params)
      cat("  maximum out-degree:                   ", x$learning$args$max.out.degree, "\n")
    if ("max.degree" %in% params)
      cat("  maximum degree:                       ", x$learning$args$max.degree, "\n")

  }#ELSE

  cat("\n")

  invisible(x)

}#PRINT.BN

