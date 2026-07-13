
bn.fit.assignment.backend = function(x, name, value) {

  # preserve the original object for subsequent sanity checks.
  to.replace = x[[name]]
  new = to.replace

  if (is(to.replace, c("bn.fit.dnode", "bn.fit.onode"))) {

    # check the consistency of the new conditional distribution.
    value = check.dnode.rvalue(value, node = name)
    # sanity check the new object by comparing it to the old one.
    value = check.rvalue.vs.dnode(value, to.replace)
    # replace the conditional probability table.
    new$prob = value

  }#THEN
  else if (is(to.replace, "bn.fit.gnode")) {

    if (is(value, c("lm", "glm", "penfit")) && is(to.replace, "bn.fit.gnode")) {

      # ordinary least squares, ridge, lasso, and elastic net.
      coef = .coefficients(value)
      resid = .residuals(value)
      fitted = .fitted(value)
      sd = cgsd(resid[!is.na(resid)], p = length(coef))

      # zap small values in low-order regressions to match fast.lm().
      if ((length(coef) <= 3) && isTRUE(all.equal(sd, 0))) {

        coef = zapsmall(coef)
        sd = 0
        resid = rep(0, length(resid))

      }#THEN

      value = list(coef = coef, resid = resid, fitted = fitted, sd = sd)
      # if the intercept is not there, set it to zero.
      if ("(Intercept)" %!in% names(value$coef))
        value$coef = c("(Intercept)" = 0, value$coef)

    }#THEN
    else {

      # check the consistency of the new conditional distribution.
      value = check.gnode.rvalue(value, node = name)

    }#ELSE

    # sanity check the new object by comparing it to the old one.
    value = check.rvalue.vs.gnode(value, to.replace)

    # replace the regression coefficients, keeping the names and the ordering.
    if (is.null(names(value$coef)))
      new$coefficients = structure(value$coef, names = names(new$coefficients))
    else
      new$coefficients = noattr(value$coef[names(new$coefficients)], ok = "names")

    # replace the residuals' standard deviation.
    new$sd = noattr(value$sd)

    # replace the residuals, padding with NAs if needed.
    if (is.null(value$resid))
      new$residuals = rep(as.numeric(NA), length(new$resid))
    else
      new$residuals = noattr(value$resid)

    # replace the fitted values, padding with NAs if needed.
    if (is.null(value$fitted))
      new$fitted.values = rep(as.numeric(NA), length(new$fitted))
    else
      new$fitted.values = noattr(value$fitted)

  }#THEN
  else if (is(to.replace, "bn.fit.cgnode")) {

    # carry discrete parents' configurations from the old object.
    value$configs = to.replace$configs
    # check the consistency of the new conditional distribution.
    value = check.cgnode.rvalue(value, node = name)
    # sanity check the new object by comparing it to the old one.
    check.rvalue.vs.cgnode(value, to.replace)

    # replace the regression coefficients, keeping the names and the ordering.
    if (is.null(names(value$coef)))
      dimnames(value$coef) = dimnames(new$coefficients)

    new$coefficients = noattr(value$coef)

    # replace the residuals' standard deviation.
    new$sd = structure(noattr(value$sd), names = colnames(value$coef))

    # replace the residuals, padding with NAs if needed.
    if (is.null(value$resid))
      new$residuals = rep(as.numeric(NA), length(new$resid))
    else
      new$residuals = noattr(value$resid)

    # replace the fitted values, padding with NAs if needed.
    if (is.null(value$fitted))
      new$fitted.values = rep(as.numeric(NA), length(new$fitted))
    else
      new$fitted.values = noattr(value$fitted)

  }#THEN
  else if (is(to.replace, "bn.fit.zihpnode")) {

    # check the consistency of the new conditional distribution.
    value = check.zihpnode.rvalue(value, node = name)
    # sanity check the new object by comparing it to the old one.
    value = check.rvalue.vs.zihpnode(value, to.replace)

    # replace the regression coefficients, keeping the names and the ordering.
    if (is.null(names(value$inflation)))
      new$inflation = structure(value$inflation, names = names(new$inflation))
    else
      new$inflation = noattr(value$inflation[names(new$inflation)], ok = "names")
    if (is.null(names(value$intensity)))
      new$intensity = structure(value$intensity, names = names(new$intensity))
    else
      new$intensity = noattr(value$intensity[names(new$intensity)], ok = "names")

    # replace the dispersion.
    new$dispersion = noattr(value$dispersion)

    # replace the residuals, padding with NAs if needed.
    if (is.null(value$resid))
      new$residuals = rep(as.numeric(NA), length(new$resid))
    else
      new$residuals = noattr(value$resid)

    # replace the fitted values, padding with NAs if needed.
    if (is.null(value$fitted))
      new$fitted.values = rep(as.numeric(NA), length(new$fitted))
    else
      new$fitted.values = noattr(value$fitted)

  }#THEN
  else if (is(to.replace, "bn.fit.zinbnode")) {

    # check the consistency of the new conditional distribution.
    value = check.zinbnode.rvalue(value, node = name)
    # sanity check the new object by comparing it to the old one.
    value = check.rvalue.vs.zinbnode(value, to.replace)

    # replace the regression coefficients, keeping the names and the ordering.
    if (is.null(names(value$inflation)))
      new$inflation = structure(value$inflation, names = names(new$inflation))
    else
      new$inflation = noattr(value$inflation[names(new$inflation)], ok = "names")
    if (is.null(names(value$prsucc)))
      new$prsucc = structure(value$prsucc, names = names(new$prsucc))
    else
      new$prsucc = noattr(value$prsucc[names(new$prsucc)], ok = "names")

    # replace the intensity.
    new$failures = noattr(value$failures)

    # replace the residuals, padding with NAs if needed.
    if (is.null(value$resid))
      new$residuals = rep(as.numeric(NA), length(new$resid))
    else
      new$residuals = noattr(value$resid)

    # replace the fitted values, padding with NAs if needed.
    if (is.null(value$fitted))
      new$fitted.values = rep(as.numeric(NA), length(new$fitted))
    else
      new$fitted.values = noattr(value$fitted)

  }#THEN
  else {

     stop("unknown node type '", class(to.replace), "'.")

  }#ELSE

  return(new)

}#BN.FIT.ASSIGNMENT.BACKEND
