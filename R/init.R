
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

  bnlearn.classes = c("bn", "bn.fit", "bn.naive", "bn.tan")

  setHook(packageEvent("graph", "attach"), action = "append",
    function(...) {

      # re-register S4 methods.
      for (cl in bnlearn.classes) {

        setMethod("nodes", cl, where = topenv(parent.frame()),
          function(object) .nodes(object))
        setMethod("nodes<-", cl, where = topenv(parent.frame()),
          function(object, value) .relabel(object, value))
        setMethod("degree", cl, where = topenv(parent.frame()),
          function(object, Nodes) .degree(object, Nodes))

      }#FOR

  })

  setHook(packageEvent("BiocGenerics", "attach"), action = "append",
    function(...) {

      # re-register S4 methods.
      for (cl in bnlearn.classes) {

        setMethod("path", cl, where = topenv(parent.frame()),
          function(object, from, to, direct = TRUE, underlying.graph = FALSE,
            debug = FALSE)
              path.exists(object, from = from, to = to, direct = direct,
                 underlying.graph = underlying.graph, debug = debug))
        setMethod("score", cl, where = topenv(parent.frame()),
          function(x, data, type = NULL, ..., by.node = FALSE, debug = FALSE)
            network.score(x = x, data = data, type = type, ...,
              by.node = by.node, debug = debug))

      }#FOR

  })

  setHook(packageEvent("igraph", "attach"), action = "append",
    function(...) {

      # re-register the S3 method for as.igraph().
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
    tryMethod("path", cl,
      definition = function(object, from, to, direct = TRUE,
                     underlying.graph = FALSE, debug = FALSE)
                       path.exists(object, from = from, to = to, direct = direct,
                          underlying.graph = underlying.graph, debug = debug),
      generic = function(object, ...) standardGeneric("path"))
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

