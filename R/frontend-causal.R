# mutilated network used for interventions and in likelihood weighting.
mutilated = function(x, evidence) {

  check.bn.or.fit(x)
  # check the evidence, disallowing non-ideal interventions if needed.
  evidence = check.evidence(evidence, graph = x,
                ideal.only = is(x, "bn.fit.gnet"))

  if (is(x, "bn")) {

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

    return(mutilated.backend.bn(x, evidence))

  }#THEN
  else {

    return(mutilated.backend.fitted(x, evidence))

  }#ELSE

}#MUTILATED

# twin network used for counterfactuals.
twin = function(x) {

  check.bn.or.fit(x)

  if (is(x, "bn.fit"))
    x = bn.net(x)
  else {

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

  }#ELSE

  twin.backend.bn(x)

}#TWIN

# perform an intervention on a network or fitted network.
intervention = function(x, evidence) {

  # the mutilation encodes the intervention into the network.
  mutilated(x = x, evidence = evidence)

}#INTERVENTION

# set up a counterfactual network.
counterfactual = function(x, evidence, merging = TRUE) {

  check.bn.or.fit(x)
  check.logical(merging)

 # check the evidence, disallowing non-ideal interventions if needed.
  evidence = check.evidence(evidence, graph = x,
                ideal.only = is(x, "bn.fit.gnet"))

  if (is(x, "bn.fit"))
    x = bn.net(x)
  else {

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

  }#ELSE

  counterfactual.backend.bn(x = x, evidence = evidence, merging = merging)

}#COUNTERFACTUAL

