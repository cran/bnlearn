
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

  # if no generic is present, create it.
  if ("graph" %!in% loadedNamespaces()) {

    setGeneric("nodes", function(object, ...) standardGeneric("nodes"))
    setGeneric("nodes<-", function(object, value) standardGeneric("nodes<-"))
    setGeneric("degree", function(object, Nodes, ...) standardGeneric("degree"))

  }#THEN

  # add the methods.
  for (cl in bnlearn.classes) {

    setMethod("nodes", cl, function(object) .nodes(object))
    setMethod("nodes<-", cl, function(object, value) .relabel(object, value))
    setMethod("degree", cl, function(object, Nodes) .degree(object, Nodes))

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
