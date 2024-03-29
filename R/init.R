
# safe version of setMethod().
tryMethod = function(f, signature, definition, generic) {

  # try a first time.
  catch = try(setMethod(f, signature, definition), silent = TRUE)

  if (is(catch, "try-error")) {

    # if it failed, create a generic function ...
    setGeneric(f, generic)
    # ... and then try again.
    setMethod(f, signature, definition)

  }#THEN

}#TRYMETHOD

# set up hooks, S4 classes and initialize global variables.
.onLoad = function(lib, pkg) {

  bnlearn.classes = c("bn", "bn.fit", available.classifiers)

  # re-register S4 methods (from .GlobalEnv, as we would do manually in an
  # interactive session.
  setHook(packageEvent("graph", "attach"), action = "append",
    function(...) {

      for (cl in bnlearn.classes) {

        setMethod("nodes", cl, where = .GlobalEnv,
          function(object) .nodes(object))
        setMethod("nodes<-", cl, where = .GlobalEnv,
          function(object, value) .relabel(object, value))
        setMethod("degree", cl, where = .GlobalEnv,
          function(object, Nodes) .degree(object, Nodes))

      }#FOR

  })

  setHook(packageEvent("gRbase", "attach"), action = "append",
    function(...) {

      for (cl in bnlearn.classes) {

        setMethod("nodes", cl, where = .GlobalEnv,
          function(object) .nodes(object))

      }#FOR

  })

  setHook(packageEvent("BiocGenerics", "attach"), action = "append",
    function(...) {

      for (cl in bnlearn.classes) {

        setMethod("score", cl, where = .GlobalEnv,
          function(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)
            network.score(x = x, data = data, type = type, ...,
              by.node = by.node, debug = debug))

      }#FOR

  })

  setHook(packageEvent("igraph", "attach"), action = "append",
    function(...) {

      registerS3method("as.igraph", class = "bn", method = as.igraph.bn,
        envir = asNamespace("igraph"))
      registerS3method("as.igraph", class = "bn.fit", method = as.igraph.bn.fit,
        envir = asNamespace("igraph"))

  })

  # make bnlearn's classes known to S4.
  setClass("bn")
  setClass("bn.fit")
  setClass("bn.naive")
  setClass("bn.tan")

  # add the methods (if no generic is present, create it) .
  for (cl in bnlearn.classes) {

    tryMethod("nodes", cl,
      definition = function(object) .nodes(object),
      generic = function(object, ...) standardGeneric("nodes"))
    tryMethod("nodes<-", cl,
      definition = function(object, value) .relabel(object, value),
      generic = function(object, value) standardGeneric("nodes<-"))
    tryMethod("degree", cl,
      definition = function(object, Nodes) .degree(object, Nodes),
      generic = function(object, Nodes, ...) standardGeneric("degree"))
    tryMethod("score", cl,
      definition = function(x, data, type = NULL, ..., by.node = FALSE,
                     debug = FALSE)
                       network.score(x = x, data = data, type = type, ...,
                         by.node = by.node, debug = debug),
      generic = function (x, ...) standardGeneric("score"))

  }#FOR

  # load the shared library.
  library.dynam("bnlearn", package = pkg, lib.loc = lib)
  # initialize stuff at the C level.
  .Call(call_onLoad)

}#.ONLOAD

# clean up global variables.
.onUnload = function(libpath) {

  # initialize stuff at the C level.
  .Call(call_onUnload)
  # unload the shared library.
  library.dynam.unload("bnlearn", libpath = libpath)

}#ON.UNLOAD

