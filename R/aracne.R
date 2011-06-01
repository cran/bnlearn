
aracne.backend = function(x, estimator, whitelist, blacklist, debug = FALSE) {

  # fix the whitelist and the blacklist to keep the C side simple.
  nodes = names(x)

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

  }#THEN

  .Call("aracne",
        data = x,
        estimator = estimator,
        whitelist = whitelist,
        blacklist = blacklist,
        debug = debug,
        PACKAGE = "bnlearn")

}#ARACNE.BACKEND

