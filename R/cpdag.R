
# reconstruct the equivalence class of a network.
cpdag.backend = function(x, moral = FALSE, fix = FALSE, wlbl = TRUE,
    debug = FALSE) {

  nodes = names(x$nodes)

  # reset wlbl if x contains no whitelist and no blacklist.
  if ((is.null(x$learning$whitelist)) && (is.null(x$learning$blacklist)))
    wlbl = FALSE

  amat = .Call(call_cpdag,
               arcs = x$arcs,
               nodes = nodes,
               moral = moral,
               fix = fix,
               wlbl = wlbl,
               whitelist = x$learning$whitelist,
               blacklist = x$learning$blacklist,
               illegal = x$learning$illegal,
               debug = debug)

  # update the arcs of the network.
  x$arcs = amat2arcs(amat, nodes)
  # update the network structure.
  x$nodes = cache.structure(nodes, amat = amat)

  return(x)

}#CPDAG.BACKEND

# backend to get a DAG out of a CPDAG (still in the same equivalence class).
cpdag.extension = function(x, debug = FALSE) {

  nodes = names(x$nodes)

  # update the arcs of the network.
  x$arcs = cpdag.arc.extension(arcs = x$arcs, nodes = nodes, debug = debug)
  # update the network structure.
  x$nodes = cache.structure(nodes, arcs = x$arcs)

  return(x)

}#CPDAG.EXTENSION

# backend to get a set of directed arcs out of a CPDAG.
cpdag.arc.extension = function(arcs, nodes, debug = FALSE) {

  .Call(call_pdag_extension,
        arcs = arcs,
        nodes = nodes,
        debug = debug)

}#CPDAG.ARC.EXTENSION

# produce all possible extensions of a CPDAG.
cextend.all.backend = function(x, debug = FALSE) {

  # if all arcs are directed, there is nothing to extend.
  nodes = names(x$nodes)
  undirected = which.undirected(x$arcs, nodes)

  if (!any(undirected))
    return(x)

  # drop all directed arcs.
  ug = empty.graph.backend(nodes)
  ug$arcs = x$arcs[undirected, , drop = FALSE]
  ug$nodes = cache.structure(nodes, arcs = ug$arcs)

  # drop all isolated nodes, which are irrelevant to the algorithm.
  isolated = (sapply(names(ug$nodes), .degree, x = ug) == 0)
  ug = subgraph.backend(ug, names(which(!isolated)))

  # reset the recursion depth limit to allow the expansion.
  recursion.depth.limit = options("expressions")
  if (recursion.depth.limit < length(ug$nodes)) {

    # the largest valid value is 500000.
    options("expressions" = min(length(ug$nodes) + 50, 500000))
    on.exit(options("expressions" = recursion.depth.limit))

  }#THEN

  # produce all possible combinations of directions of the undirected arcs.
  extensions = uccg.all.extensions(ug, debug = debug)

  # put back the directed arcs.
  extensions = lapply(extensions, function(ext) {

    dag = empty.graph.backend(nodes)
    dag$arcs = rbind(ext, x$arcs[!undirected, , drop = FALSE])
    dag$nodes = cache.structure(nodes, arcs = dag$arcs)

    return(dag)

  })

  return(extensions)

}#CEXTEND.ALL.BACKEND

# enumerate all possible extensions of an undrected connected chordal graph.
uccg.all.extensions = function (ug, state, generate = TRUE, debug = FALSE) {

  n = nnodes(ug)

  # initialise the global state at the top of the recursion.
  if (missing(state)) {

    state = new.env()
    state$A = vector(n, mode = "list")
    state$A[[1]] = nodes(ug)

    state$tau = character(0)

    state$count = 0
    state$AMOs = vector(0, mode = "list")

  }#THEN

  # tau contains a complete ordering, produce the directed extension.
  if (length(state$tau) == n) {

    state$count = state$count + 1

    if (debug)
      cat("  @ found consistent extension #", state$count, "!\n")

    if (generate) {

      state$AMOs =
        c(state$AMOs, list(pdag2dag.backend(ug$arcs, ordering = state$tau)))

    }#THEN

    return(NULL)

  }#THEN

  # here starts the actual recursive exploration.
  i = max(which(sapply(state$A, length) > 0))
  v = state$A[[i]][1]
  x = v

  repeat {

    # add node x to the node ordering, and remove it from the candidates.
    state$A[[i]] = setdiff(state$A[[i]], x)
    state$tau = c(state$tau, x)

    if (debug)
      cat("* considering ordering:", state$tau, "\n")

    for (w in setdiff(ug$nodes[[x]]$nbr, state$tau)) {

      j = which(sapply(state$A, `%in%`, x = w))
      state$A[[j]] = setdiff(state$A[[j]], w)
      state$A[[j + 1]] = c(state$A[[j + 1]], w)

    }#FOR

    # recurse.
    uccg.all.extensions(ug, state, generate = generate, debug = debug)

    # remove node x from the node ordering.
    for (w in setdiff(ug$nodes[[x]]$nbr, state$tau)) {

      j = which(sapply(state$A, `%in%`, x = w))
      state$A[[j]] = setdiff(state$A[[j]], w)
      state$A[[j - 1]] = c(state$A[[j - 1]], w)

    }#FOR

    state$A[[i]] = c(state$A[[i]], x)
    state$tau = setdiff(state$tau, x)

    # update the node search space.
    if (x == v) {

      subg = subgraph.backend(ug, state$A[[i]])
      R = sapply(setdiff(names(subg$nodes), v), path.exists, x = subg, from = v)

      if (length(R) == 0)
        R = character(0)
      else
        R = names(which(R))

    }#THEN

    if (length(R) == 0)
      break
    else {

      x = R[length(R)]
      R = R[-length(R)]

    }#ELSE

  }#REPEAT

  if (generate)
    return(state$AMOs)
  else
    return(state$count)

}#UCCG.ALL.EXTENSIONS
