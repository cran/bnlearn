
# backend for conditional probability queries.
conditional.probability.query = function(fitted, event, evidence, method,
    extra, debug = FALSE) {

  if (method == "ls") {

    logic.sampling(fitted = fitted, event = event, evidence = evidence,
      n = extra$n, batch = extra$batch, debug = debug)

  }#THEN

}#CONDITIONAL.PROBABILITY.QUERY

logic.sampling = function(fitted, event, evidence, n, batch, debug) {

  cpxe = cpe = 0L
  filtered = logical(n)
  matching = logical(n)
  r = logical(n)

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (n in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and 
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (n > 0) {

      if (is.fitted.discrete(fitted))
        generated.data = rbn.discrete(x = fitted, n = n)
      else
        generated.data = rbn.continuous(x = fitted, n = n)

    }#THEN
    else
      break
  
    if (debug)
      cat("* generated", n, "samples from the bayesian network.\n")
  
    # evaluate the expression defining the evidence.
    r = eval(evidence, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != n)
      stop("logical vector for evidence is of length ", length(r), 
        " instead of ", n, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)
 
    # update the global counters.
    cpe = cpe + length(which(filtered))
 
    if (debug) {

      lwfilter = length(which(filtered))
      cat("  > evidence matches", lwfilter, "samples out of", n,
        "(p =", lwfilter/n, ").\n");

    }#THEN

    # evaluate the expression defining the event.
    r = eval(event, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != n)
      stop("logical vector for event is of length ", length(r), 
        " instead of ", n, ".")
    # filter out the samples not matching the event we are looking for.
    matching = filtered & r & !is.na(r)
 
    # update the global counters.
    cpxe = cpxe + length(which(matching))
 
    if (debug) {

      lwmatch = length(which(matching))
      lwratio = ifelse(lwfilter == 0, 0, lwmatch/lwfilter)
      cat("  > event matches", lwmatch, "samples out of", lwfilter, "(p =",
        lwratio, ").\n");

    }#THEN

  }#FOR

  # prevent divide-by-zero errors.
  result = ifelse(cpe == 0, 0, cpxe / cpe)

  if (debug && (nbatches > 1)) {

    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event matches", cpxe, "samples out of", cpe, "(p =", result, ").\n");

  }#THEN

  return(result)

}#LOGIC.SAMPLING

