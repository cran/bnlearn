.onLoad <- function(lib, pkg) {

    # silence all warnings while looking for suggested packages.
    warning.level  = as.numeric(options("warn"))
    options("warn" = -1)

    # load graohviz and set the corresponding flag.
    graphviz.enabled <<- require(Rgraphviz)

    if (!graphviz.enabled)
      if (!("Rgraphviz" %in% rownames(installed.packages())))
        cat("Package Rgraphviz not installed.\n")
      else
        cat("Package Rgraphviz not loaded.\n")

    # restore the original warning level.
    options("warn" = warning.level)

    # load the shared library.
    library.dynam("bnlearn", package = pkg, lib.loc = lib)

}

