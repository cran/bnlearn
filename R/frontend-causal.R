
# convert bayesian networks into structural causal models.
as.scm = function(x, ...) {

  UseMethod("as.scm")

}#AS.SCM

as.scm.bn = function(x, ...) {

  check.bn(x)

  # the network must be a directed acyclic graph to treat it as causal.
  if (!is.completely.directed(x))
    stop("the graph is only partially directed.")
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the graph contains cycles.")

  from.bn.to.scm(x)

}#AS.SCM.BN

# convert bayesian networks into structural models.
as.bn.scm = function(x, ...) {

  check.scm(x)

  from.scm.to.bn(x)

}#AS.BN.SCM

# twin network used for counterfactuals.
twin = function(x) {

  if (missing(x) || !is(x, c("bn", "scm", "bn.fit")))
    stop("an object of class 'bn', 'scm' or 'bn.fit' is required.")

  if (is(x, "scm")) {

    twin = twin.backend(x)

  }#THEN
  else if (is(x, "bn")) {

    # already a twin network, nothing to do.
    if (is(x, "bn.twin"))
      return(x)

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

    twin = twin.backend(from.bn.to.scm(x))
    twin = from.scm.to.bn(twin)

  }#THEN
  else if (is(x, "bn.fit")) {

    if (is(x, "bn.fit.gnet"))
      twin = twin.backend.fitted(x)
    else {

      # fall back: return only the twin network for other bn.fit objects.
      scm = from.bn.to.scm(bn.net(x))
      twin = twin.backend(scm)
      twin = from.scm.to.bn(twin)

    }#ELSE

  }#THEN

  return(twin)

}#TWIN

# perform an intervention on a network or fitted network.
intervention = function(x, evidence) {

  if (missing(x) || !is(x, c("bn", "scm", "bn.fit")))
    stop("an object of class 'bn', 'scm' or 'bn.fit' is required.")
  # check the evidence, disallowing non-ideal interventions if needed.
  evidence = check.evidence(evidence, graph = x,
               ideal.only = is(x, c("bn.fit.gnet", "bn.fit.cgnet")))

  if (is(x, "scm")) {

    scm = intervention.backend(x, evidence = evidence)

  }#THEN
  else if (is(x, "bn")) {

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

    scm = intervention.backend(from.bn.to.scm(x), evidence = evidence)
    scm = from.scm.to.bn(scm)

  }#THEN
  else if (is(x, "bn.fit")) {

    scm = intervention.backend.fitted(x, evidence = evidence)

  }#THEN

  return(scm)

}#INTERVENTION

# mutilated() is an alias of intervention().
mutilated = intervention

# set up a counterfactual network.
counterfactual = function(x, evidence, merging = TRUE) {

  if (missing(x) || !is(x, c("bn", "scm", "bn.fit")))
    stop("an object of class 'bn', 'scm' or 'bn.fit' is required.")
  if (is(x, c("bn.ctf", "scm.ctf")))
    stop("'x' is already a counterfactual network.")
  check.logical(merging)

  if (is(x, "scm")) {

    twin = twin.backend(x)

    # check the evidence, disallowing non-ideal interventions if needed.
    evidence = check.evidence(evidence, graph = twin, ideal.only = FALSE,
                 counterfactual = TRUE)

    ctf = counterfactual.backend(twin, evidence = evidence, merging = merging)

  }#THEN
  else if (is(x, "bn")){

    # the network must be a directed acyclic graph to treat it as causal.
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")
    if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
      stop("the graph contains cycles.")

    scm = from.bn.to.scm(x)
    twin = twin.backend(scm)

    # check the evidence, disallowing non-ideal interventions if needed.
    evidence = check.evidence(evidence, graph = twin, ideal.only = FALSE,
                 counterfactual = TRUE)

    ctf = counterfactual.backend(twin, evidence = evidence, merging = merging)
    ctf = from.scm.to.bn(ctf)

  }#ELSE
  else if (is(x, "bn.fit")) {

    if (is(x, "bn.fit.gnet")) {

      twin = twin.backend.fitted(x)

      # check the evidence, disallowing non-ideal interventions if needed.
      evidence = check.evidence(evidence, graph = twin, ideal.only = TRUE,
                   counterfactual = TRUE)

      ctf = counterfactual.backend.fitted(twin, evidence = evidence,
              merging = merging)

    }#THEN
    else {

      scm = from.bn.to.scm(bn.net(x))
      twin = twin.backend(scm)

      # check the evidence, disallowing non-ideal interventions if needed.
      evidence = check.evidence(evidence, graph = twin, ideal.only = FALSE,
                   counterfactual = TRUE)

      ctf = counterfactual.backend(twin, evidence = evidence, merging = merging)
      ctf = from.scm.to.bn(ctf)

    }#ELSE

  }#THEN

  return(ctf)

}#COUNTERFACTUAL

