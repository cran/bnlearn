.onLoad <- function(lib, pkg) {

    # silence all warnings while looking for suggested packages.
    warning.level  = as.numeric(options("warn"))
    options("warn" = -1)

    # require the utils package explicitly.
    require(utils)
    # load graohviz and set the corresponding flag.
    if ("Rgraphviz" %in% rownames(installed.packages()))
      graphviz.enabled <<- require(Rgraphviz)

    cat("Package Rgraphviz", 
      ifelse(graphviz.enabled, "loaded successfully.\n", "not loaded.\n"))

    # restore the original warning level.
    options("warn" = warning.level)

    # load the shared library.
    library.dynam("bnlearn", package = pkg, lib.loc = lib)

}#.ONLOAD

