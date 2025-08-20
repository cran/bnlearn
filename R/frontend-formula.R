
# describe the network with a "model string".
modelstring = function(x) {

  check.bn.or.fit(x)
  # no model string if the graph is partially directed.
  if (is(x, "bn"))
    if (!is.completely.directed(x))
      stop("the graph is only partially directed.")

  modelstring.backend(x)

}#MODELSTRING

# set a specific network structure with the model string.
"modelstring<-" = function(x, debug = FALSE, value) {

  check.modelstring(value)
  check.logical(debug)

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
model2network = function(string, ordering = NULL, debug = FALSE) {

  # check string's class and format.
  check.modelstring(string)
  # check the node ordering; NULL is ok this time, it lets the backend decide.
  if (!is.null(ordering))
    check.nodes(ordering)

  check.logical(debug)

  result = model2network.backend(string, node.order = ordering, debug = debug)

  # check the node ordering again now that the graph is built.
  if (!is.null(ordering)) {

    check.nodes(ordering, graph = result, min.nodes = length(result$nodes),
      max.nodes = length(result$nodes))

  }#THEN

  return(result)

}#MODEL2NETWORK

# model-string-to-bn conversion function.
as.bn.character = function(x, ...) {

  # warn about unused arguments.
  check.unused.args(list(...), character(0))

  model2network(x)

}#AS.BN.CHARACTER

