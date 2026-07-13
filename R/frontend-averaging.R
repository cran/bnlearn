
# build the averaged network structure using arc strengths and a
# significance threshold.
averaged.network = function(strength, threshold) {

  # check the main argument.
  check.bn.strength(strength, valid = c("bootstrap", "bayes-factor"))
  # check the strength threshold.
  threshold = check.threshold(threshold, strength = strength)

  avg = averaged.network.backend(strength = strength, threshold = threshold)

  # add the metadata for the print() method.
  avg$learning$algo = "averaged"
  avg$learning$args = list(threshold = threshold)

  return(avg)

}#AVERAGED.NETWORK

# compute the inclusion threshold from a set of arc strengths.
inclusion.threshold = function(strength) {

  # check the (only) argument.
  check.bn.strength(strength, valid = c("bootstrap", "bayes-factor"))

  threshold(strength = strength, method = "l1")

}#INCLUSION.THRESHOLD

# average (the parameters of) multiple bn.fit objects with identical structures.
mean.bn.fit = function(x, ..., weights = NULL) {

  # check the bn.fit objects.
  fitted = c(list(x), list(...))
  # check the nodes are the same, and that the object has the right structure.
  for (s in seq_along(fitted))
    check.fitted.vs.fitted(fitted[[1]], fitted[[s]])
  # check the weights.
  weights = check.weights(weights, length(fitted))

  # average the objects.
  average.fitted(fitted, weights)

}#MEAN.BN.FIT

