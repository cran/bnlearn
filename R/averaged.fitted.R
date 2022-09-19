
# average discrete nodes whose parameters are organised in a conditional
# probability table.
mean.fitted.dnode = function(node, fitted, weights) {

  # create a conditional probability table filled with zeroes...
  cpt = fitted[[1]][[node]]$prob
  cpt[] = 0
  # ... add up the conditional probability tables from the networks, with
  # the respective weights...
  for (i in seq_along(fitted))
    cpt = cpt + fitted[[i]][[node]]$prob * weights[i]
  # ... and normalize the result so that columns sum up to 1 again.
  cpt = cpt / sum(weights)

}#MEAN.FITTED.DNODE

# average continuous nodes whose parameters are organised in vectors or matrices
# of regression coefficients and standard errors.
mean.fitted.cgnode = function(node, fitted, weights) {

  # allocate the vector of the regression coefficients and the standard
  # error, both zeroed...
  coefs = fitted[[1]][[node]]$coefficients
  coefs[] = 0
  sd = fitted[[1]][[node]]$sd
  sd[] = 0
  # ... add up both of them, with the respective weights...
  for (i in seq_along(fitted)) {

    coefs = coefs + fitted[[i]][[node]]$coefficients * weights[i]
    sd = sd + fitted[[i]][[node]]$sd * weights[i]

  }#FOR
  # ... and normalize the result.
  coefs = coefs / sum(weights)
  sd = sd / sum(weights)

  return(list(coef = coefs, sd = sd))

}#MEAN.FITTED.CGNODE

# average several bn.fit objects with the same structure.
mean.fitted = function(fitted, weights) {

  # all the networks have the same structure and parameter sets, so we can
  # allocate the return value by copying one of them.
  averaged = fitted[[1]]
  cl = class(averaged)
  class(averaged) = "list"

  for (node in names(averaged)) {

    if (is(averaged[[node]], c("bn.fit.dnode", "bn.fit.onode"))) {

      averaged[[node]]$prob =
        mean.fitted.dnode(node = node, fitted = fitted, weights = weights)

    }#THEN
    else if (is(averaged[[node]], c("bn.fit.gnode", "bn.fit.cgnode"))) {

      averaged[[node]][c("coefficients", "sd")] =
        mean.fitted.cgnode(node = node, fitted = fitted, weights = weights)

      # in addition to averaging the parameters, remove the fitted values and
      # the residuals if present.
      averaged[[node]]$fitted.values = as.numeric(NA)
      averaged[[node]]$residuals = as.numeric(NA)
      # same with the configurations of the discrete parents in a conditional
      # Gaussian node.
      if ("configs" %in% names(averaged[[node]]))
        averaged[[node]]$configs =
          factor(NA, levels = levels(averaged[[node]]$configs))

    }#THEN

  }#FOR

  class(averaged) = cl

  # preserve attributes for classifiers when all fitted networks have the same
  # type and training variable.
  if (all(sapply(fitted, is, available.classifiers))) {

    classifier.type =
      sapply(fitted, function(x) intersect(class(x), available.classifiers))
    classifier.training = sapply(fitted, attr, "training")

    if ((length(unique(classifier.type)) == 1) &&
        (length(unique(classifier.training)) == 1)) {

      class(averaged) = union(unique(classifier.type), class(averaged))
      attr(averaged, "training") = unique(classifier.training)

    }#THEN

  }#THEN

  return(averaged)

}#MEAN.FITTED
