
# describe the network with a "model string".
modelstring = function(x) {

  # check x's class.
  check.bn(x)
  # no model string if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  formula.backend(x)

}#MODELSTRING

"modelstring<-" = function(x, value) {

  res = model2network(value)

  # check that the new network constains the same nodes as the old one.
  if (!setequal(names(x$nodes), names(res$nodes)))
    stop("the model string describes a different network (the nodes are different).")

  return(res)

}#MODELSTRING<-

# bn-to-character (i.e. the model string) conversion function.
# an alias of modelstring().
as.character.bn = function(x, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  modelstring(x)

}#AS.CHARACTER.BN

# generate an object of class bn from a model string.
model2network = function(string, debug = FALSE) {

  # check string's class.
  if (!is(string, "character"))
    stop("string must be a character string.")
  # check debug.
  check.logical(debug)

  model2network.backend(string, debug = debug)

}#MODEL2NETWORK

# model-string-to-bn conversion function.
as.bn.character = function(x, debug = FALSE) {

  model2network(x, debug = debug)

}#AS.BN.CHARACTER

