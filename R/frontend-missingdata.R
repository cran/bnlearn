
# impute missing data from a bn.fit object.
impute = function(object, data, method, ..., debug = FALSE) {

  # check the data are there.
  check.data(data, allow.levels = TRUE, allow.missing = TRUE,
    warn.if.no.missing = TRUE)
  # check whether the data agree with the bayesian network.
  check.fit.vs.data(object, data)
  # check the imputation method.
  method = check.imputation.method(method, data)
  # check the extra arguments passed to the imputation methods.
  extra.args = check.imputation.extra.args(method, list(...))
  # check debug.
  check.logical(debug)

  impute.backend(fitted = object, data = data, method, extra.args = extra.args,
    debug = debug)

}#IMPUTE

# structural expectation-maximization.
structural.em = function(x, maximize = "hc", maximize.args = list(), fit = "mle",
    fit.args = list(), impute, impute.args = list(), return.all = FALSE,
    start = NULL, max.iter = 5, debug = FALSE) {

  ntests = 0

  # check the data are there.
  data.info =
    check.data(x, allow.levels = TRUE, allow.missing = TRUE,
      warn.if.no.missing = TRUE, stop.if.all.missing = !is(start, "bn.fit"))

  # check the max.iter parameter.
  max.iter = check.max.iter(max.iter)
  # check debug and return.data.
  check.logical(debug)
  check.logical(return.all)

  # check the arguments used for structure learning as in rsmax2().
  check.learning.algorithm(algorithm = maximize, class = "score")

  critical.arguments = c("x", "heuristic", "start", "debug")
  check.unused.args(intersect(critical.arguments, names(maximize.args)),
    character(0))
  maximize.args[critical.arguments] =
    list(x = NULL, heuristic = maximize, start = NULL, debug = debug)

  # check the arguments used for parameter learning as in bn.cv().
  check.fitting.method(method = fit, data = x)
  fit.args = check.fitting.args(method = fit, network = NULL, data = x,
               extra.args = fit.args)

  # check the arguments used for imputations.
  impute = check.imputation.method(impute, x)
  impute.args = check.imputation.extra.args(impute, impute.args)

  # if there is no preseeded network, use an empty one.
  if (is.null(start)) {

    dag = empty.graph(nodes = names(x))
    fitted = bn.fit.backend(dag, data = x, method = fit, extra.args = fit.args,
               data.info = data.info)

  }#THEN
  else {

    # check start's class.
    check.bn.or.fit(start)
    # check the preseeded network against the data set.
    if (is(start, "bn")) {

      check.bn.vs.data(start, x)
      # check the preseeded network is valid for the model assumptions.
      check.arcs.against.assumptions(start$arcs, x, maximize.args$score)

      dag = start
      fitted = bn.fit.backend(start, data = x, method = fit,
                 extra.args = fit.args, data.info = data.info)

    }#THEN
    else if (is(start, "bn.fit")) {

      check.fit.vs.data(start, x)
      # fitted networks are necessarily valid for the model assumptions.

      dag = bn.net(start)
      fitted = start

    }#THEN

  }#ELSE

  # initialize the algorithm.
  if (debug) {

    cat("* initializing the network to perform the first imputation.\n")
    cat("* network structure:\n")
    print(dag)
    cat("* fitted parameters:\n")
    print(fitted)

  }#THEN

  # from this point on the data are complete.
  data.info$complete.nodes[names(x)] = TRUE

  for (i in seq(max.iter)) {

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* iteration", i, ", expectation step .\n")

    }#THEN

    # expectation step.
    complete = impute.backend(fitted = fitted, data = x, method = impute,
                 extra.args = impute.args, debug = debug)

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* iteration", i, ", maximization step .\n")

    }#THEN

    # maximization step, structure learning (starting from the previous network).
    maximize.args$x = complete
    maximize.args$start = dag
    dag = do.call("greedy.search", maximize.args)

    # maximization step, parameter learning.
    fitted.new = bn.fit.backend(dag, data = complete, method = fit,
                   extra.args = fit.args, data.info = data.info)

    if (debug) {

      cat("----------------------------------------------------------------\n")
      cat("* fitted parameters:\n")
      print(fitted)

    }#THEN

    # update the number of model comparisons.
    ntests = ntests + dag$learning$ntests

    # a rudimentary stopping rule (the fitted network must be the same as in
    # the previous iteration).
    if (isTRUE(all.equal(fitted, fitted.new)))
      break
    else
      fitted = fitted.new

  }#FOR

  # set the metadata.
  dag$learning$algo = "sem"
  dag$learning$maximize = maximize
  dag$learning$impute = impute
  dag$learning$fit = fit
  dag$learning$ntests = ntests

  if (return.all)
    invisible(list(dag = dag, imputed = complete, fitted = fitted))
  else
    invisible(dag)

}#STRUCTURAL.EM

