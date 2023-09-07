
# sanitize the extra arguments passed to the random graph generation algorithms.
check.graph.generation.args = function(method, nodes, extra.args) {

  # check the probability of generating arcs in a random graph.
  if (has.argument(method, "prob", graph.generation.extra.args))
    extra.args[["prob"]] =
      check.arc.generation.prob(extra.args[["prob"]], nodes = nodes)

  # check the thinning factor of a Markov chain.
  if (has.argument(method, "every", graph.generation.extra.args))
    extra.args[["every"]] =
      check.thinning(extra.args[["every"]])

  # this magic number for the burn-in comes from the reference implementation
  # of Ide-Cozman.
  if (has.argument(method, "burn.in", graph.generation.extra.args))
    extra.args[["burn.in"]] =
      check.burn.in(extra.args[["burn.in"]], default = 6 * length(nodes)^2)

  if (has.argument(method, "max.in.degree", graph.generation.extra.args))
    extra.args[["max.in.degree"]] =
      check.max.degree(extra.args[["max.in.degree"]], label = "in-degree",
        nodes = nodes)

  if (has.argument(method, "max.out.degree", graph.generation.extra.args))
    extra.args[["max.out.degree"]] =
      check.max.degree(extra.args[["max.out.degree"]], label = "out-degree",
        nodes = nodes)

  if (has.argument(method, "max.degree", graph.generation.extra.args)) {

    extra.args[["max.degree"]] =
      check.max.degree(extra.args[["max.degree"]], label = "degree",
        nodes = nodes)

    if (is.finite(extra.args[["max.in.degree"]]) &&
        extra.args[["max.in.degree"]] > extra.args[["max.degree"]])
      stop("the maximun in-degree must be lesser or equal to the maximum degree.")

    if (is.finite(extra.args[["max.out.degree"]]) &&
        extra.args[["max.out.degree"]] > extra.args[["max.degree"]])
      stop("the maximun out-degree must be lesser or equal to the maximum degree.")

  }#THEN

  # warn about and remove unused arguments.
  extra.args = check.unused.args(extra.args, graph.generation.extra.args[[method]])

  return(extra.args)

}#CHECK.GRAPH.GENERATION.ARGS

# check the probability of generating arcs in a random graph.
check.arc.generation.prob = function(prob, nodes) {

  if (!is.null(prob)) {

    # prob must be a probability.
    if (!is.probability(prob))
      stop("the branching probability must be a numeric value in [0,1].")

  }#THEN
  else {

    # this default produces graphs with about the same number of arcs as there
    # are nodes.
    prob = 2 / (length(nodes) - 1)

  }#ELSE

  return(prob)

}#CHECK.ARC.GENERATION.PROB

# check the thinning factor of a Markov chain.
check.thinning = function(every) {

  if (!is.null(every)) {

    if (!is.positive(every))
      stop("'every' must be a positive integer number.")

  }#THEN
  else {

    # default to no thinning.
    every = 1

  }#ELSE

  return(every)

}#CHECK.THINNING

# check the length of burn-in of a Markov chain.
check.burn.in = function(burn.in, default) {

  if (!is.null(burn.in)) {

    if (!is.positive(burn.in))
      stop("the burn in length must be a positive integer number.")

  }#THEN
  else {

    burn.in = default

  }#ELSE

  return(burn.in)

}#CHECK.BURN.IN

# check bounds on various types of node degrees.
check.max.degree = function(degree, label, nodes) {

  if (!is.null(degree)) {

    if (!is.positive.integer(degree))
      stop("the maximum ", label, " must be a positive integer number.")

    if (degree >= length(nodes)) {

      warning("a node cannot have ", label,
              " greater or equal to the number of nodes in the graph.")
      warning("the condition on the ", label, " will be ignored.")

    }#THEN

  }#THEN
  else {

    degree = Inf

  }#ELSE

  return(degree)

}#CHECK.MAX.DEGREE

