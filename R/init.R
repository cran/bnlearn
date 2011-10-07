.onLoad = function(lib, pkg) {

    # silence all warnings while looking for suggested packages.
    warning.level  = as.numeric(options("warn"))
    options("warn" = -1)

    graphviz.enabled = lattice.enabled = FALSE

    # require the utils package explicitly.
    require(utils)
    # load graohviz and set the corresponding flag.
    if ("Rgraphviz" %in% rownames(installed.packages())) {

      require(grid)
      require(graph)

      graphviz.enabled = TRUE

    }#THEN

    packageStartupMessage("Package Rgraphviz ",
      ifelse(graphviz.enabled, "will be loaded as needed.", "not loaded."))

    # load lattice and set the corresponding flag.
    if ("lattice" %in% rownames(installed.packages()))
      lattice.enabled = require(lattice)

    packageStartupMessage("Package lattice ",
      ifelse(lattice.enabled, "loaded successfully.", "not loaded."))

    # restore the original warning level.
    options("warn" = warning.level)

    # load the shared library.
    library.dynam("bnlearn", package = pkg, lib.loc = lib)

}#.ONLOAD

