
# describe the network with a "model string".
modelstring = function(x) {

  # check x's class.
  check.bn(x)
  # no model string if the graph is partially directed.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")

  formula.backend(x)

}#MODELSTRING

# set a specific network structure with the model string.
"modelstring<-" = function(x, debug = FALSE, value) {

  # check string's class.
  if (!is(value, "character"))
    stop("the model string must be a character string.")

  model2network.backend(value, node.order = names(x$nodes), debug = debug)

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

# generic method to fool R CMD check (it's not even exported).
as.bn = function(string, ...) {
 
   UseMethod(string)
 
}#AS.BN

# model-string-to-bn conversion function.
as.bn.character = function(string, debug = FALSE) {

  model2network(string, debug = debug)

}#AS.BN.CHARACTER

