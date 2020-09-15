
# perform conditional probability queries.
cpquery = function(fitted, event, evidence, cluster = NULL, method = "ls", ...,
    debug = FALSE) {

  # check fitted's class.
  check.fit(fitted)
  # check debug.
  check.logical(debug)
  # check event and evidence.
  if (missing(event))
    stop("the expression describing the event is missing.")
  if (missing(evidence))
    stop("the expression describing the evidence is missing.")
  # check the generation method.
  check.label(method, choices = cpq.algorithms, labels = cpq.labels,
    argname = "query method", see = "cpquery")
  # check the cluster.
  if (!is.null(cluster)) {

    check.cluster(cluster)

    # disable debugging, the slaves do not cat() here.
    if (debug) {

      warning("disabling debugging output for parallel computing.")
      debug = FALSE

    }#THEN

  }#THEN

  extra.args = check.cpq.args(fitted = fitted, event = NULL,
                 extra.args = list(...), method = method, action = "cpquery")

  # deparse the expression for the event before passing it to
  # the backend and beyond.
  event = substitute(event)

  # recheck event and evidence expression after deparsing.
  if (!(is.language(event) || identical(event, TRUE)))
    stop("event must be an unevaluated expression or TRUE.")
  if (method == "ls") {

    if (missing(evidence))
      stop("the expression describing the evidence is missing.")

    # deparse evidence expression before passing it to the backend and beyond.
    evidence = substitute(evidence)
    # recheck event and evidence expression after deparsing.
    if (!(is.language(evidence) || identical(evidence, TRUE)))
      stop("evidence must be an unevaluated expression or TRUE.")

  }#THEN
  else if (method == "lw") {

    evidence = check.evidence(evidence, fitted)

  }#THEN

  # special-case pointless queries.
  if (isTRUE(event))
    return(1)

  conditional.probability.query(fitted = fitted, event = event,
    evidence = evidence, method = method, extra = extra.args, cluster = cluster,
    debug = debug)

}#CPQUERY

