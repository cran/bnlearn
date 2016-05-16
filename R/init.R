
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

  }#FOR

  # load the shared library.
  library.dynam("bnlearn", package = pkg, lib.loc = lib)
  # initialize stuff at the C level.
  .Call("c_onLoad")

}#.ONLOAD

# clean up global variables.
.onUnload = function(libpath) {

  # initialize stuff at the C level.
  .Call("c_onUnload")
  # unload the shared library.
  library.dynam.unload("bnlearn", libpath = libpath)

}#ON.UNLOAD
