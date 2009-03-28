
random.graph.backend = function(num, nodes, method, extra.args, debug) {

  if (method == "ordered") {

    res = ordered.graph(num = num, nodes = nodes, prob = extra.args$prob)

  }#THEN
  else if (method == "ic-dag") {

    res = ide.cozman.graph(num = num, nodes = nodes,
            burn.in = extra.args$burn.in,
            max.in.degree = extra.args$max.in.degree,
            max.out.degree = extra.args$max.out.degree,
            max.degree = extra.args$max.degree,
            debug = debug)

  }#THEN
  else if (method == "empty") {

    res = empty.graph.backend(num = num, nodes = nodes)

  }#THEN

  return(res)

}#RANDOM.GRAPH.BACKEND

# generate a random directed acyclic graph.
ordered.graph = function (num, nodes, prob) {

  .Call("ordered_graph",
        nodes = nodes,
        num = as.integer(num),
        prob = prob,
        PACKAGE = "bnlearn")

}#ORDERED.GRAPH

# generate a random directed acyclic graph accordin to a uniform
# probability distribution over the space of connected graphs.
ide.cozman.graph = function(num, nodes, burn.in, max.in.degree,
    max.out.degree, max.degree, debug) {

  .Call("ide_cozman_graph",
        nodes = nodes,
        num = as.integer(num),
        burn.in = as.integer(burn.in),
        max.in.degree = as.numeric(max.in.degree),
        max.out.degree = as.numeric(max.out.degree),
        max.degree = as.numeric(max.degree),
        debug = debug,
        PACKAGE = "bnlearn")

}#IDE.COZMAN.GRAPH

# generate an empty 'bn' object given a set of nodes.
empty.graph.backend = function(nodes, num = 1) {

  .Call("empty_graph", 
        nodes = nodes,
        num = as.integer(num),
        PACKAGE = "bnlearn")

}#EMPTY.GRAPH.BACKEND
