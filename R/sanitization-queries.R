
# sanitize the extra arguments passed to the conditional probability algorithms.
check.cpq.args = function(fitted, event, extra.args, method, action) {

  if (method %in% c("ls", "lw")) {

    if (!is.null(extra.args$n)) {

      if (!is.positive.integer(extra.args$n))
        stop("the number of observations to be sampled must be a positive integer number.")

    }#THEN
    else {

      # this is a rule of thumb, the error of the estimate has no closed-form
      # expression (Koller & Friedman).
      if (!is(fitted, "bn.fit.gnet"))
        extra.args$n = 5000 * max(1, round(log10(nparams.fitted(fitted))))
      else
        extra.args$n = 500 * nparams.fitted(fitted)

    }#ELSE

    if (!is.null(extra.args$batch)) {

      if ((action == "cpdist") && (method == "lw")) {

        extra.args$batch = NULL
        warning(" 'batch' will be ignored for speed and memory efficience.")

      }#THEN
      else {

        if (!is.positive.integer(extra.args$batch))
          stop("the number of observations to be sampled must be a positive integer number.")

        if (extra.args$batch > extra.args$n) {

          warning("cannot generate a batch bigger than the whole generated data set.")
          warning("batch size will be ignored.")

        }#THEN

      }#ELSE

    }#THEN
    else {

      # perform small simulations in a single batch, and split larger ones.
      extra.args$batch = min(extra.args$n, 10^4)

    }#ELSE

    if (!is.null(extra.args$query.nodes)) {

      check.nodes(extra.args$query.nodes, graph = fitted)

      # make sure the nodes to be simulated are included.
      if (action == "cpdist")
        extra.args$query.nodes = c(event, extra.args$query.nodes)

    }#THEN

  }#THEN

  # warn about unused arguments.
  check.unused.args(extra.args, cpq.extra.args[[method]])

  return(extra.args)

}#CHECK.CPQ.ARGS

# check evidence in list format for mutilated networks.
check.mutilated.evidence = function(evidence, graph) {

  # check whether evidence is there.
  if (missing(evidence))
    stop("evidence must be a list with elements named after the nodes in the graph.")
  # if evindence is TRUE there's nothing to check.
  if (identical(evidence, TRUE))
    return(TRUE)
  # check whether evidence is a named list.
  if (!is(evidence, "list"))
    stop("evidence must be a list with elements named after the nodes in the graph.")
  # check the node labels in evidence.
  check.nodes(names(evidence), graph = graph)

  if (is(graph, "bn")) {

    # check the network is completely directed.
    if (!is.dag(graph$arcs, names(graph$nodes)))
      stop("the graph is only partially directed.")

  }#THEN
  else {

     # check the evidence is appropriate for the nodes.
     for (fixed in names(evidence)) {

       # extract the node and the evidence.
       cur = graph[[fixed]]
       ev = evidence[[fixed]]

       if (is(cur, c("bn.fit.dnode", "bn.fit.onode"))) {

         if (is.factor(ev))
           evidence[[fixed]] = ev = as.character(ev)

         if (!is.string.vector(ev) || any(ev %!in% dimnames(cur$prob)[[1]]))
           stop("the evidence for node ", fixed, " must be valid levels.")

       }#THEN
       else if (is(cur, "bn.fit.gnode")) {

         # for continuous nodes evidence must be real numbers.
         if (!is.real.vector(ev) || (length(ev) %!in% 1:2))
           stop("the evidence ", fixed, " must be a real number or a finite interval.")
         storage.mode(ev) = "double"
         # make sure interval boundaries are in the right order.
         evidence[[fixed]] = sort(ev)

       }#THEN

     }#FOR

  }#THEN

  return(evidence)

}#CHECK.MUTILATED.EVIDENCE

