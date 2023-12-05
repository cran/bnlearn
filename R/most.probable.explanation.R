
# partition a query into independent subqueries using the d-separation induced
# by the conditioning on the event nodes.
query.partitioning = function(fitted, event, evidence, debug = FALSE) {

  catchall = list(event = character(0), evidence = evidence)
  queries = list()

  dag = bn.net(fitted)

  # until all the event nodes in the query have been evaluated...
  while (length(event) > 0) {

    # ... take the first one...
    target = event[1]
    # ... find out its markov blanket...
    target.mb = dag$nodes[[target]]$mb

    if (debug)
      cat("  > considering event node", target, "with markov blanket",
          target.mb, "\n")

    # ... while the markov blanket(s) of the event node(s) contain other
    # event nodes, merge them...
    while (length(intersect(event, target.mb)) > 0) {

      if (debug)
        cat("    > nodes", intersect(event, target.mb),
            "are also event nodes, considering them as well.\n")

      target = c(target, intersect(event, target.mb))
      target.mb = unique(unlist(sapply(target, function(t) dag$nodes[[t]]$mb)))
      target.mb = setdiff(target.mb, target)

      if (debug) {

        cat("    > the joint markov blanket of", target, "is",
            ifelse(length(target.mb) == 0, "empty",
                   paste(target.mb, collapse = " ")), ".\n")

      }#THEN

    }#THEN

    if (all(target.mb %in% evidence)) {

      # ... remove the node from consideration...
      event = setdiff(event, target)
      # ... and store the subquery.
      queries = c(queries, list(list(event = target, evidence = target.mb)))

    }#THEN
    else {

      # ... remove the node from consideration...
      event = setdiff(event, target)
      # ... and give up because the markov blanket is not observable because of
      # the limited set of evidence nodes.
      catchall$event = c(catchall$event, target)

    }#ELSE

  }#WHILE

  # in the simple case where the catchall has a single event node, try to reduce
  # it using d-separation.

  # merge back events for which reduction to a subquery was impossible...
  if (length(catchall$event) > 0) {

    # ... but try to reduce the evidence first using d-separation if there is a
    # single event node left.
    if (length(catchall$event) == 1) {

      for (e in catchall$evidence)
        if (dseparation(dag, catchall$event, e, setdiff(catchall$evidence, e)))
          catchall$evidence = setdiff(catchall$evidence, e)

    }#THEN

    queries = c(queries, list(catchall))

  }#THEN

  if (debug) {

    cat("  @ query partitioned into:\n")
    for (q in queries) {

      cat("    event:", q$event, "\n    evidence:",
          ifelse(length(q$evidence) == 0, "(empty)",
                 paste(q$evidence, collapse = " ")), "\n")

    }#FOR

  }#DEBUG

  return(queries)

}#QUERY.PARTITIONING

# most probable explanation for discrete networks, with exact inference.
mpe.discrete.query = function(jtree, event, evidence, value) {

  # if the network is not a junction tree yet, transform it...
  if (is(jtree, "bn.fit"))
    jtree = from.bn.fit.to.grain(jtree, compile = TRUE)

  # ... incorporate the evidence, if any...
  if (length(evidence) == 0) {

    jpred = jtree

  }#THEN
  else {

    jpred = gRain::setEvidence(jtree, nodes = evidence,
              states = sapply(value[, evidence, drop = FALSE], as.character))

    # ... give up if it is not possible to observe the evidence...
    if (gRain::pEvidence(jpred) <= sqrt(.Machine$double.eps))
      return(NULL)

  }#ELSE

  # ... get the joint probability table of the event nodes...
  ppp = grain.query(jpred, nodes = event, type = "joint")
  # ... find the coordinates of the largest probability...
  id = which(ppp == max(ppp), arr.ind = TRUE)
  # ... pick one at random if there are multiple maxima...
  id = id[sample(nrow(id), size = 1), , drop = FALSE]
  # ... and replace them with the corresponding levels.
  mpe = id[, event, drop = FALSE]
  for (node in colnames(mpe))
   mpe[1, node] = dimnames(ppp)[[node]][id[1, node]]

  return(mpe)

}#MPE.DISCRETE.QUERY

# most probable explanation for gaussian networks, with exact inference.
mpe.gaussian.query = function(mvn, event, evidence, value) {

  if (length(evidence) == 0) {

    # if there is no evidence, use the marginal expectations of the event
    # variables.
    mpe = mvn$mu[event]

  }#THEN
  else {

    # if there is evidence, use the exectation of the conditional distribution
    # of the event nodes given the evidence nodes.
    mpe = conditional.mvnorm(mu = mvn$mu, sigma = mvn$sigma, to = event,
            from = evidence, value = value[, evidence])

  }#ELSE

  return(mpe)

}#MPE.GAUSSIAN.QUERY

