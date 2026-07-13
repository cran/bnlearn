
chow.liu.backend = function(data, nodes, estimator, whitelist, blacklist,
    conditional = NULL, debug = FALSE) {

  # fix the whitelist and the blacklist to keep the C side simple.
  if (!is.null(blacklist)) {

    # arcs must be blacklisted in both directions, so keep only
    # the undirected ones.
    blacklist = blacklist[which.undirected(blacklist, nodes), , drop = TRUE]
    # keep only one direction for each blacklisted arc.
    blacklist = pdag2dag.backend(blacklist, nodes)

  }#THEN

  if (!is.null(whitelist)) {

    # keep only one direction for each whitelisted arc.
    whitelist = pdag2dag.backend(whitelist, nodes)

    # the chow-liu algorithms allows the selection of exactly length(nodes) arcs,
    # so the whitelist must contain less.
    if (nrow(whitelist) > length(nodes))
      stop("too many whitelisted arcs, there can be only ", length(nodes), ".")

  }#THEN

  .Call(call_chow_liu,
        data = data,
        nodes = nodes,
        estimator = estimator,
        whitelist = whitelist,
        blacklist = blacklist,
        complete = attr(data, "metadata")$complete.nodes,
        conditional = conditional,
        debug = debug)

}#CHOW.LIU.BACKEND

