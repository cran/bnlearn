
# load suggested packages and initialize global variables.
.onLoad = function(lib, pkg) {

  # set the test/score counter.
  assign(".test.counter", 0, envir = .GlobalEnv)

  # load the shared library.
  library.dynam("bnlearn", package = pkg, lib.loc = lib)

}#.ONLOAD

