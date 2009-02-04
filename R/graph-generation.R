
random.graph.backend = function(num, nodes, method, extra.args, debug) {

  if (method == "ordered") {

    res = ordered.graph(num = num, nodes = nodes, prob = extra.args$prob, 
            debug = debug)

  }#THEN
  else if (method == "ic-dag") {

    res = ide.cozman.graph(num = num, nodes = nodes, 
            burn.in = extra.args$burn.in, 
            max.in.degree = extra.args$max.in.degree,
            max.out.degree = extra.args$max.out.degree,
            max.degree = extra.args$max.degree,
            debug = debug)

  }#THEN
  else if (method == "empty") {

    em = empty.graph.backend(nodes)
    res = sapply(1:num, function(x) { em }, simplify = FALSE)

  }#THEN

  if (num == 1)
    return(res[[1]])
  else
    return(res)

}#RANDOM.GRAPH.BACKEND

# generate a random directed acyclic graph.
ordered.graph = function (num, nodes, prob, debug) {

    i = 1
    n = length(nodes)
    amat = matrix(0L, ncol = n, nrow = n, dimnames = list(nodes, nodes))
    mask = upper.tri(amat)

    res = list()
    dummy = empty.graph.backend(nodes)
    dummy$learning$algo = "ordered"
    dummy$learning$args$prob = prob

    repeat {

      # sample each arc with probability 'prob'.
      amat[mask] = sample(c(0L, 1L), n * (n - 1) / 2, 
                     prob = c(1 - prob, prob), replace = TRUE)
      # build the corresponding object of class bn.
      dummy$arcs = amat2arcs(amat, nodes)
      dummy$nodes = cache.structure(nodes, amat = amat)

      if (debug)
        cat(formula.backend(dummy), "\n")

      res[[i]] = dummy

      if (i < num)
        i = i + 1
      else
        break

    }#REPEAT

return(res)

}#ORDERED.GRAPH

ide.cozman.graph = function(num, nodes, burn.in, max.in.degree, 
    max.out.degree, max.degree, debug) {

  i = 1
  n = length(nodes)
  to.be.updated = TRUE

  # initialize the result list.
  res = list()
  # initialize the dummy bn object with all the metadata.
  dummy = empty.graph.backend(nodes)
  dummy$learning$algo = "ic-dag"
  dummy$learning$args$burn.in = burn.in
  dummy$learning$args$max.in.degree = max.in.degree
  dummy$learning$args$max.out.degree = max.out.degree
  dummy$learning$args$max.degree = max.degree

  # initialize a simple ordered tree with n nodes, where all nodes 
  # have just one parent, except the first one that does not have 
  # any parent.
  amat = matrix(0L, ncol = n, nrow = n)
  amat[cbind(1:(n-1),  2:n)] = 1L

  repeat {

    arc = sample(n, 2, replace = FALSE)    

    if (debug) {

      cat("* current model (", i, ") is:\n")
      print(amat2arcs(amat, nodes))

    }#THEN

    if (amat[arc[1], arc[2]] > 0) {

      # if the arc (i, j) exists in the actual graph, delete the arc, 
      # provided that the underlying graph remains connected.

      if (debug)
        cat("  > arc", nodes[arc[1]], "->", nodes[arc[2]], "is present.\n")

      # if there is a(n undirected) path in the underlying (undirected) graph
      # from arc[1] to arc[2], the graph is still connected.
      if (has.path(arc[1], arc[2], 1:n, amat, underlying.graph = TRUE, 
            exclude.direct = TRUE)) {

        if (debug)
          cat("  @ removing arc", nodes[arc[1]], "->", nodes[arc[2]], ".\n")

        amat[arc[1], arc[2]] = 0L
        to.be.updated = TRUE

      }#THEN
      else {

        if (debug)
          cat("  > not removing arc", nodes[arc[1]], "->", nodes[arc[2]], " (graph not connected!).\n")

      }#ELSE

    }#THEN
    else {

      # add the arc, provided that the underlying graph remains acyclic.

      if (debug)
        cat("  > arc", nodes[arc[1]], "->", nodes[arc[2]], "is not present.\n")

      # if there is a (directed) path from arc[2] to arc[1], adding arc[1] -> arc[2]
      # would create a cycle, so do not do it.
      if (!has.path(arc[2], arc[1], 1:n, amat)) {

        # check the conditions on the {in-,out-,}degrees of the nodes.
        if ((any(in.degree(amat) >= max.in.degree)) ||
            (any(out.degree(amat) >= max.out.degree)) ||
            (any(degree(amat) >= max.degree))) {

          if (debug)
            cat("  > not adding arc", nodes[arc[1]], "->", nodes[arc[2]], " (degree!).\n")

        }#THEN
        else {

          if (debug)
            cat("  @ adding arc", nodes[arc[1]], "->", nodes[arc[2]], ".\n")

          amat[arc[1], arc[2]] = 1L
          to.be.updated = TRUE

        }#ELSE

      }#THEN
      else {

        if (debug)
          cat("  > not adding arc", nodes[arc[1]], "->", nodes[arc[2]], " (cycles!).\n")

      }#THEN

    }#ELSE

    if (i < burn.in) {

      # if this iteration is still in the "burn in", increment the counter
      # but throw away the current model.
      i = i + 1

    }#THEN
    else if (i >= burn.in && i < burn.in + num) {

      # update the network structure.
      if (to.be.updated) {

        # save the new arc set.
        dummy$arcs = amat2arcs(amat, nodes)

        # check which nodes have to be updated.
        if (i == burn.in)
          changed.nodes = nodes
        else
          changed.nodes = unique(c(nodes[arc], dummy$nodes[[arc[1]]]$mb, dummy$nodes[[arc[2]]]$mb))

        if (debug)
          cat("  > updating nodes", changed.nodes, ".\n")

        # update the chosen nodes.
        for (updated in changed.nodes)
          dummy$nodes[[updated]] = cache.partial.structure(nodes, target = updated, amat = amat)
        # reset the to.be.updated flag.
        to.be.updated = FALSE

      }#THEN

      # the markov chain is now stationary, so this model is a good one;
      # add it to the list to be returned.
      res[[i - burn.in + 1]] = dummy
      i = i + 1

    }#THEN
    else
      break

  }#REPEAT

  return(res)

}#IDE.COZMAN.GRAPH

# generate an empty 'bn' object given a set of nodes.
empty.graph.backend = function(nodes) {

  arcs = matrix(character(0), ncol = 2, dimnames = list(c(), c("from", "to")))
  learning.structure = structure(lapply(nodes, function(n) {
        list(mb = character(0), nbr = character(0)) }), names = nodes)
  nodes.structure = structure(lapply(nodes, function(n) {
        list(mb = character(0), nbr = character(0), parents = character(0),
          children = character(0)) }), names = nodes)
  res = structure(list(
    learning = list(
      nodes = learning.structure,
      arcs = arcs,
      whitelist = NULL,
      blacklist = NULL,
      test = "none",
      ntests = 0,
      algo = "empty",
      args = list()
    ),
    nodes = nodes.structure,
    arcs = arcs
  ), class = "bn")

  res

}#EMPTY.GRAPH.BACKEND
