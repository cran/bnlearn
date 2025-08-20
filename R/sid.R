# compute the structural interventional distance between two DAGs.
sid.dag.vs.dag = function(learned, true, debug = FALSE) {

  incorrect = 0
  nodes = nodes(true)
  pmat = path.matrix(true)

  # iterate over the intervention target nodes...
  for (target in nodes) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* checking causal effects from node", target, ".\n")

    }#THEN

    true.target.parents = parents(true, target)
    learned.target.parents = parents(learned, target)
    true.target.children = children(true, target)

    reachable = causal.reachability(dag = true, pmat = pmat, target = target,
                  adjustment.set = learned.target.parents, debug = debug)

    # ... and all the other nodes to check the interventional distributions.
    for (other in setdiff(nodes, target)) {

      ijGNull = ijHNull = FALSE

      # the causal effect is zero in the true model.
      if (pmat[target, other] == FALSE)
        ijGNull = TRUE

      # the causal effect is zero in the learned model, because the intervention
      # removed the target node parents when creating the mutilated graph.
      if (other %in% learned.target.parents)
        ijHNull = TRUE

      # if the causal effect is zero in the learned network but not in the true
      # one, the effect of the intervention does not match.
      if (ijHNull && !ijGNull) {

        if (debug)
          cat("@ interventional distributions are different for node", other, ".\n")

        incorrect = incorrect + 1
        next

      }#THEN

      # if both conditions hold, the effect of the intervention matches.
      if ((ijHNull && ijGNull) ||
          setequal(learned.target.parents, true.target.parents) ) {

        next

      }#THEN

      # children of the target node in the true network that have the other node
      # as a descendant.
      on.causal.path = names(which(pmat[true.target.children, other]))

      if (any(pmat[on.causal.path, learned.target.parents])) {

        if (debug)
          cat("@ interventional distributions are different for node", other, ".\n")

        incorrect = incorrect + 1

      }#THEN

      # checking paths that are not completely directed from the targe to the
      # other node but are nevertheless open because of v-structures.
      if (reachable[match(other, nodes)]) {

        if (debug)
          cat("@ interventional distributions are different for node", other, ".\n")

        incorrect = incorrect + 1

      }#THEN

    }#FOR

  }#FOR

  return(incorrect)

}#SID.DAG.VS.DAG

# similar to an adjacency matrix, but encodes directed paths instead.
path.matrix = function(dag, debug = FALSE) {

  .Call(call_path_matrix,
        x = dag,
        debug = debug)

}#PATH.MATRIX

# similar to an adjacency matrix, but encodes d-separations instead.
causal.reachability = function(dag, pmat, target, adjustment.set, debug = FALSE) {

  p = nnodes(dag)
  top = 1:p
  bottom = p + 1:p
  reach.pa = reach.ch = structure(rep(FALSE, p), names = nodes(dag))

  # parents are reachable non-causal paths.
  reach.pa[parents(dag, target)] = TRUE

  rmat = .Call(call_reachability_matrix,
               x = dag,
               path.matrix = pmat,
               target.node = target,
               adjustment.set = adjustment.set,
               debug = debug)

  # other nodes reachable through the parents or children of the target node.
  paths = colSums(rmat[bottom, ][reach.pa, , drop = FALSE])
  reach.pa[names(which(paths > 0))] = TRUE
  paths = colSums(rmat[top, ][reach.ch, , drop = FALSE])
  reach.ch[names(which(paths > 0))] = TRUE

  return(reach.ch | reach.pa)

}#CAUSAL.REACHABILITY

