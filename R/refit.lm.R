
# backend for the as.lm() functions.
lm.refit.node = function(node, data) {

  # construct the model formula, inserting the intercept term as needed.
  if (length(node$parents) == 0)
    model = paste(node$node, "~ 1")
  else
    model = paste(node$node, "~", paste(node$parents, collapse = "+"))
  # fit the model.
  lm.fit = lm(model, data = data)
  # replace the formula in the call to lm() for legibility.
  lm.fit$call$formula = formula(model)

  return(lm.fit)

}#LM.REFIT.NODE

