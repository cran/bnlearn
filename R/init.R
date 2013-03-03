
# load suggested packages and initialize global variables.
.onLoad = function(lib, pkg) {

  setHook(packageEvent("graph", "attach"), action = "append",
    function(...) {

      # re-register S4 methods.
      setMethod("nodes", "bn", where = topenv(parent.frame()),  
        function(object) .nodes(object))
      setMethod("nodes", "bn.fit", where = topenv(parent.frame()),
        function(object) .nodes(object))
      setMethod("nodes", "bn.naive", where = topenv(parent.frame()),
        function(object) .nodes(object))
      setMethod("nodes", "bn.tan", where = topenv(parent.frame()),
        function(object) .nodes(object))

      setMethod("degree", "bn", where = topenv(parent.frame()), 
        function(object, Nodes) .degree(object, Nodes))
      setMethod("degree", "bn.fit", where = topenv(parent.frame()), 
        function(object, Nodes) .degree(object, Nodes))
      setMethod("degree", "bn.naive", where = topenv(parent.frame()), 
        function(object, Nodes) .degree(object, Nodes))
      setMethod("degree", "bn.tan", where = topenv(parent.frame()), 
        function(object, Nodes) .degree(object, Nodes))

  })

  # load the shared library.
  library.dynam("bnlearn", package = pkg, lib.loc = lib)

}#.ONLOAD

