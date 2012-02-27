
random.graph.backend = function(num, nodes, method, extra.args, debug = FALSE) {

  if (method == "ordered") {

    res = ordered.graph(num = num, nodes = nodes, prob = extra.args$prob)

  }#THEN
  else if (method == "ic-dag") {

    # adjust the number of graph to generate with the stepping factor.
    num = num * extra.args$every

    res = ide.cozman.graph(num = num, nodes = nodes,
            burn.in = extra.args$burn.in,
            max.in.degree = extra.args$max.in.degree,
            max.out.degree = extra.args$max.out.degree,
            max.degree = extra.args$max.degree,
            connected = TRUE, debug = debug)

    # keep only every k-th network.
    if (num > 1) {

      res = res[seq(from = extra.args$every, to = num, by = extra.args$every)]

    }#THEN

  }#THEN
  else if (method == "melancon") {

    # adjust the number of graph to generate with the stepping factor.
    num = num * extra.args$every

    res = ide.cozman.graph(num = num, nodes = nodes,
            burn.in = extra.args$burn.in,
            max.in.degree = extra.args$max.in.degree,
            max.out.degree = extra.args$max.out.degree,
            max.degree = extra.args$max.degree,
            connected = FALSE, debug = debug)

    # keep only every k-th network.
    if (num > 1) {

      res = res[seq(from = extra.args$every, to = num, by = extra.args$every)]

    }#THEN

  }#THEN
  else if (method == "empty") {

    res = empty.graph.backend(num = num, nodes = nodes)

  }#THEN

  return(res)

}#RANDOM.GRAPH.BACKEND

# generate a random directed acyclic graph.
ordered.graph = function(num, nodes, prob) {

  .Call("ordered_graph",
        nodes = nodes,
        num = as.integer(num),
        prob = prob,
        PACKAGE = "bnlearn")

}#ORDERED.GRAPH

# generate a random directed acyclic graph accordin to a uniform
# probability distribution over the space of connected graphs (if
# connected = TRUE) or the space of graphs (if connected = FALSE).
ide.cozman.graph = function(num, nodes, burn.in, max.in.degree,
    max.out.degree, max.degree, connected, debug = FALSE) {

  .Call("ide_cozman_graph",
        nodes = nodes,
        num = as.integer(num),
        burn.in = as.integer(burn.in),
        max.in.degree = as.numeric(max.in.degree),
        max.out.degree = as.numeric(max.out.degree),
        max.degree = as.numeric(max.degree),
        connected = connected,
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
