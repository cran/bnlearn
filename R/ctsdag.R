
# fix an arc if any arc pointing to its tail is fixed, to avoid introducing
# new v-structures.
determine.arc.status = function(arc.status, cur) {

  upstream = which(arc.status[, "to"] == arc.status[cur, "from"])

  if (length(upstream) > 0) {

    upstream.arcs = arc.status[upstream, , drop = FALSE]

    if (any(upstream.arcs$fixed))
      arc.status[cur, "fixed"] = TRUE

  }#THEN

  return(arc.status)

}#DETERMINE.ARC.STATUS

# if an arc has no incoming incident arcs, apply determine.arc.status, else
# traverse the network recursively until reaching an arc with no incoming
# incident arcs that are untraversed.
arc.traverser = function(arc.status, cur) {

  # this arc has already been traversed, nothing to do.
  if (arc.status$traversed[cur])
    return(arc.status)

  upstream = which(arc.status[, "to"] == arc.status[cur, "from"])

  if (length(upstream) > 0) {

    untraversed = !arc.status[upstream, "traversed"]

    if (any(untraversed)) {

      untraversed.arcs = upstream[untraversed]

      for (u in untraversed.arcs)
        arc.status = arc.traverser(arc.status, u)

    }#THEN

    arc.status = determine.arc.status(arc.status, cur)

  }#THEN
  else {

    arc.status[cur, "traversed"] = TRUE

  }#ELSE

  return(arc.status)

}#ARC.TRAVERSER

# get any arc probabilities set by graphical priors.
get.arc.prior.probabilities = function(learning.args, nodes) {

  if (is.null(learning.args$prior) || is.null(learning.args$beta)) {

    # no per-arc prior probabilities, return an empty data frame.
    probs = data.frame(from = character(0), to = character(0),
              prob = numeric(0), stringsAsFactors = FALSE)

  }#THEN
  else if (learning.args$prior == "marginal") {

    inclusion.prob = noattr(learning.args$beta)
    all.arcs = expand.grid(to = nodes, from = nodes, stringsAsFactors = FALSE)
    all.arcs = all.arcs[all.arcs$from != all.arcs$to, ]
    all.arcs = data.frame(all.arcs, prob = inclusion.prob)

    probs = data.frame(from = character(0), to = character(0),
              prob = numeric(0), stringsAsFactors = FALSE)

  }#THEN
  else if (learning.args$prior == "vsp") {

    warning("the 'vsp' prior does not decompose into individual arc probabilities, ignoring.")

    probs = data.frame(from = character(0), to = character(0),
              prob = numeric(0), stringsAsFactors = FALSE)

  }#THEN
  else if (learning.args$prior == "cs") {

    # the Castelo & Siebes prior is already saved in the right format.
    probs = learning.args$beta

  }#THEN

  if (!is(probs, "prior.cs"))
    probs = cs.completed.prior(probs, nodes, learning = FALSE)

  return(probs)

}#GET.ARC.PRIOR.PROBABILITIES

# Identify which edges in beta argument have non-equivalent
# prior orientation probabilities.
# Returns:
#   (logical): Logical vector correspond to arcs in beta, where element is TRUE if arc
#       has unequal direction probabilities, FALSE otherwise.
skewed.beta = function(beta, dag, debug = FALSE) {

  default_prob_arcs = !which.listed(beta[beta$prob == 1/3, , drop = F], arcs(dag), both = TRUE)

  if (any(default_prob_arcs) & debug) {

    cat("----------------------------------------------------------------\n")
    cat(" At least one arc in beta arg has prob of 1/3 and no specified reverse prob.\n
        CS prior will assign default reverse prob of 1/3 -- i.e. PDAG algorithm will assume directions equi-probable.\n")

  }#THEN

  node_names = nodes(dag)
  beta_cs = cs.completed.prior(beta, node_names)
  near_prob = apply(beta_cs[, c("fwd", "bkwd")], 1, function(row) {
    (abs(row["fwd"] - row["bkwd"]) < 0.01) &
    (round(row["fwd"], 7) != round(row["bkwd"], 7))
  })

  if (any(near_prob) & debug) {

    cat("----------------------------------------------------------------\n")
    cat(" One or more arcs in beta argument have nearly equal probability of being in either direction, consider making them equal. See ?ctsdag.\n")

  }#THEN

  matching_prob_cs_dex = apply(beta_cs[, c("fwd", "bkwd")], 1, function(row) {
    round(row["fwd"], 7) == round(row["bkwd"], 7) # 7th decimal captures irrational number equivalence
  })

  matching_prob_cs_arcs = beta_cs[matching_prob_cs_dex, c("from", "to"), drop = FALSE]
  matching_prob_dex = which.listed(as.matrix(beta[, c("from", "to")]),
                                   as.matrix(matching_prob_cs_arcs),
                                   either = TRUE)
  skewed_prob_dex = !matching_prob_dex

  return(skewed_prob_dex)

}#SKEWED.BETA

# Modify beta with whitelist and blacklist arguments.
modify.beta.wlbl = function(dag, beta, wl = NULL, bl = NULL, debug = FALSE) {

  beta.arcs = as.matrix(beta[, c("fwd", "bkwd")])

  if (!is.null(wl)) {

    wl = as.matrix(wl)
    directed = which.directed(wl)
    directed.wl = wl[directed, , drop = FALSE]
    undirected.wl = wl[!directed, , drop = FALSE]

    # arcs that are whitelisted in a single direction have all prior probability
    # assigned to that direction.
    match = which.listed(beta.arcs, directed.wl)
    beta[match, c("fwd", "bkwd")] = c(1, 0)
    match = which.listed(beta.arcs, directed.wl[, c("to", "from"), drop = FALSE])
    beta[match, c("fwd", "bkwd")] = c(0, 1)

    # arcs that are whitelisted in both directions have prior probability 0.5
    # for each direction.
    match = which.listed(beta.arcs, undirected.wl, either = TRUE)
    beta[match, c("fwd", "bkwd")] = c(0.5, 0.5)

  }#THEN

  if (!is.null(bl)) {

    bl = as.matrix(bl)
    directed = which.directed(bl)
    directed.bl = bl[directed, , drop = FALSE]
    undirected.bl = bl[!directed, , drop = FALSE]

    # disregard arcs that are blacklisted in one direction because they are
    # whitelisted in the other direction.
    if (!is.null(wl)) {

      other = which.listed(directed.bl, directed.wl[, c("to", "from"), drop = FALSE])
      directed.bl = directed.bl[!other, , drop = FALSE]

    }#THEN

    # arcs that are blacklisted in one direction have prior probability 0.5
    # in the other direction.
    match = which.listed(beta.arcs, directed.bl)
    beta[match, c("fwd", "bkwd")] = c(0, 0.5)
    match = which.listed(beta.arcs, directed.bl[, c("to", "from"), drop = FALSE])
    beta[match, c("fwd", "bkwd")] = c(0.5, 0)
    # arcs that are blacklisted in both direction have zero prior probability
    # for either direction.
    match = which.listed(beta.arcs, undirected.bl, either = TRUE)
    beta[match, c("fwd", "bkwd")] = c(0, 0)

  }#THEN

  return(beta)

}#MODIFY.BETA.WLBL

# find the transition sequence equivalent class.
ctsdag = function(x, exp, learning = FALSE, debug = FALSE) {

  # check the network.
  check.bn(x)
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  if (!is.acyclic(x$arcs, names(x$nodes), directed = TRUE))
    stop("the graph contains cycles.")
  # check the nodes that are the targets of the interventions.
  if (missing(exp))
    exp = character(0)
  else
    check.nodes(exp, graph = x)
  # check learning and debug.
  check.logical(learning)
  check.logical(debug)

  if (learning) {

    # take into consideration the interventions, whitelist and blacklist used
    # to learn the network structure.
    info = x$learning
    interventions = intersect(c(exp, names(info$args$exp)), nodes)

    # also take into consideration graphical priors, if any, and check whether
    # directions have the same probability or not.
    beta = get.arc.prior.probabilities(info$args, nodes)

    # adjust arc probabilities following whitelists and blacklists, which
    # effectively set some of those probabilities to zero or one.
    if (!is.null(info$whitelist) || !is.null(info$blacklist))
      beta = modify.beta.wlbl(x, beta, info$whitelist, info$blacklist, debug)

browser()

  }#THEN
  else {

    interventions = intersect(exp, nodes)
    beta = NULL

  }#ELSE

  arc.status = data.frame(
    arcs,
    skewed_prior = rep(FALSE, narcs), # direction probs are not equal. includes v-structures
    from_targeted = rep(FALSE, narcs), # from node targeted by intervention
    to_targeted = rep(FALSE, narcs), # to node targeted by intervention
    fixed = rep(FALSE, narcs), # fixed arcs will not be rendered undirected
    traversed = rep(FALSE, narcs), # used in checking to avoid introducing new v-structures
    stringsAsFactors = FALSE
  )

  # Construct an edge list of arcs belonging to immoral v-stuctures
  vstruct_arcs = colliders.backend(x, return.arcs = TRUE,
                   including.unshielded = TRUE, including.shielded = FALSE)

  if (nrow(vstruct_arcs) > 0)
    arc.status[which.listed(arcs, vstruct_arcs), "skewed_prior"] = TRUE

  arc.status$from_targeted = arc.status$from %in% interventions
  arc.status$to_targeted = arc.status$to %in% interventions

  if (!is.null(beta)) {

    beta_id = skewed.beta(beta, x, debug)

    if (any(beta_id)) {

      skewed_arcs = as.matrix(beta[beta_id, -3])
      arc.status[which.listed(arcs, skewed_arcs), ]$skewed_prior = TRUE

    }#THEN

  }#THEN

  arc.status$fixed = as.logical(
    arc.status$from_targeted +
      arc.status$to_targeted +
      arc.status$skewed_prior
  )

  # Traversal determines which arcs to fix.
  # So mark arcs that are already fixed as traversed.
  arc.status$traversed = arc.status$fixed

  # Propagate orientation
  for (i in seq(nrow(arc.status)))
    arc.status = arc.traverser(arc.status, i)

  tsdag = x
  not_fixed_arcs = arc.status[!arc.status$fixed, , drop = FALSE]

  if (nrow(not_fixed_arcs) > 0) {

    for (i in 1:nrow(not_fixed_arcs)) {

      #set.edge makes edges undirected
      tsdag = set.edge(tsdag, not_fixed_arcs$from[i], not_fixed_arcs$to[i])

    }#THEN

  }#THEN

  tsdag$constraints = list(exp = interventions, beta = beta)

  return(tsdag)

}#CTSDAG
