# do a partial ordering of the nodes of a graph.
topological.ordering = function(x, start = NULL, reverse = FALSE,
    debug = FALSE) {

  if (is.null(start))
    roots = root.leaf.nodes(x, leaf = reverse)
  else
    roots = start

  to.do = .Call(call_topological_ordering,
                bn = x,
                root.nodes = roots,
                reverse = reverse,
                debug = debug)

  if (is.null(start))
    return(names(sort(to.do)))
  else
    return(names(sort(to.do[to.do > 0])))

}#SCHEDULE

# use the Logic Sampling (LS) algorithm as described in "Bayesian Artificial
# Intelligence", Korb & Nicholson, chap 3.6.1.
rbn.backend = function(x, n, fix = TRUE, debug = FALSE) {

  .Call(call_rbn_master,
        fitted = x,
        n = as.integer(n),
        fix = fix,
        debug = debug)

}#RBN.BACKEND

