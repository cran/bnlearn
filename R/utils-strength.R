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
strength2lwd = function(strength, threshold, cutpoints, method, arcs = NULL,
    debug = FALSE) {

  if (debug) {

    cat("* using threshold:", threshold, "\n")
    cat("* reported arc strength are:\n")
    if (!is.null(arcs))
      print(data.frame(arcs, strength))
    else
      print(strength)

  }#THEN

  if (method == "test") {

    lwds = lwds.from.pvalues(threshold = threshold, pvalues = strength,
             cutpoints = cutpoints)

  }#ELSE
  if (method %in% c("bootstrap", "bayes-factor")) {

    # probabilities work like p-values, only reversed.
    lwds = lwds.from.pvalues(threshold = 1 - threshold, pvalues = 1 - strength,
             cutpoints = cutpoints)

  }#THEN
  else if (method == "score") {

    lwds = lwds.from.scores(threshold = threshold, deltas = strength,
             cutpoints = cutpoints)

  }#THEN

  # arcs beyond the significance threshold are given a negative weight,
  # so that graphviz.backend() will draw them as dashed lines.
  lwds$values[lwds$values == 1] = -1

  if (debug) {

    cat("* using cut points for strength intervals:\n")
    print(lwds$cutpoints)
    cat("* arc line widths:", lwds$values, "\n")

  }#THEN

  return(lwds$values)

}#STRENGTH2LWD

# line widths from p-values: significant ones are below the threshold.
lwds.from.pvalues = function(threshold, pvalues, cutpoints) {

  # cover domain boundaries as special cases:
  #   1) if threshold == 0, only p-values exactly equal to zero are significant.
  #   2) if threshold == 1, all p-values are significant.
  if (threshold == 0)
    return(list(values = 1 + 5 * (pvalues == 0), cutpoints = 0))
  if (threshold == 1)
    return(list(values = rep(6, length(pvalues)), cutpoints = 1))

  # p-values are already scaled, so the raw quantiles are good cut points.
  lwds = cut(pvalues, cutpoints, labels = FALSE, include.lowest = TRUE)
  # reverse the levels from cut, since lower is stronger.
  lwds = length(cutpoints) - lwds

  return(list(values = lwds, cutpoints = cutpoints))

}#LWDS.FROM.PVALUES

# line widths from score deltas: significant ones are below the threshold.
lwds.from.scores = function(threshold, deltas, cutpoints) {

  # cover +/-Inf as special cases:
  #   1) if threshold == -(Inf), all deltas are significant.
  #   2) if threshold == -(-Inf), no finite deltas are significant.
  if (threshold == Inf)
    return(list(values = rep(6, length(deltas)), cutpoints = c(-Inf, Inf)))
  if (threshold == -Inf)
    return(list(values = 1 + 5 * (deltas == -Inf), cutpoints = c(-Inf, Inf)))

  lwds = cut(deltas, cutpoints, labels = FALSE)
  # reverse the levels from cut, since lower is stronger.
  lwds = length(cutpoints) - lwds

  return(list(values = lwds, cutpoints = cutpoints))

}#LWDS.FROM.SCORES
