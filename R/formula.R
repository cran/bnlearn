
# write the model formula of an object of class 'bn'.
# (ported from the deal package)
formula.backend = function(x) {

  paste(sapply(schedule(x), 
    function(node) { 

      paste("[", node, ifelse(length(x$nodes[[node]]$parents) > 0, "|", ""), 
        paste(x$nodes[[node]]$parents, sep="", collapse=":"), "]", sep = "")  

    } ), collapse = "")

}#FORMULA.BACKEND

# create an object of class 'bn' from a model formula.
# (ported from the deal package)
model2network.backend = function(modelstring, debug = FALSE) {

    nodes = c()

    if (debug)
      cat("* processing model string:\n ", modelstring, "\n")

    # split the model strings into the single nodes' strings.
    st <- strsplit(strsplit(modelstring,"\\[")[[1]],"\\]")

    arcs = lapply(st, function(xp) {

      if (debug)
       cat("  > processing node string:", xp, "\n")

      # separate the node from its parents.
      pa =  strsplit(xp[1], "\\|")[[1]]

      assign('nodes', c(nodes, pa[1]), envir = sys.frame(-2))

      # if there are parents ...
      if (length(pa) > 1) {

        prnts = strsplit(pa[2], ":")[[1]]

        if (debug)
          cat("    > detected parents '", prnts, "' for node", pa[1], "\n")

        return(cbind(from = prnts, to = pa[1]))

      }#THEN

    })

    arcs = do.call(rbind, arcs)

    res = list(learning = list(nodes = NULL, arcs = arcs, whitelist = NULL, 
            blacklist = NULL, test = "mi", alpha = 0.05, ntests = 0, 
            algo = "gs"), nodes = cache.structure(sort(nodes), arcs, 
            debug = debug), arcs = arcs)

    # fill res$learning$nodes using the relevant fields from res$nodes.
    for (node in sort(nodes)) {

      res$learning$nodes[[node]]$mb = res$nodes[[node]]$mb
      res$learning$nodes[[node]]$nbr = res$nodes[[node]]$nbr

    }#THEN

    class(res) = "bn"

    res

}#MODEL2NETWORK

