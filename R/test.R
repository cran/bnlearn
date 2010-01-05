
conditional.test = function(x, y, sx, data, test, B, learning = TRUE) {

  if (learning) {

    # update the test counter.
    assign(".test.counter", get(".test.counter", envir = .GlobalEnv) + 1,
      envir = .GlobalEnv)

  }#THEN

  sx = sx[sx != ""]
  ndata = nrow(data)
  df = NULL

  if (length(sx) == 0) {

    # Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi") {

      statistic = mi.test(data[,x], data[,y], ndata, gsquare = TRUE)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Akaike Information Criterion-like test (binary, no p-value!)
    else if (test == "aict") {

      statistic = mi.test(data[,x], data[,y], ndata, gsquare = FALSE) <
               (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) / ndata
      p.value = as.integer(statistic)

    }#THEN
    # Shrinked Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi-sh") {

      statistic = shmi.test(data[,x], data[,y], ndata, gsquare = TRUE)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (chi-square asymptotic distribution)
    else if (test == "x2") {

      statistic = x2.test(data[,x], data[,y], ndata)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Canonical (Linear) Correlation (Student's t distribution)
    else if (test == "cor") {

      statistic = fast.cor(data[,x], data[,y], ndata)
      df = ndata - 2
      p.value = pt(abs((statistic * sqrt(ndata - 2) / sqrt(1 - statistic^2))),
                  df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

       statistic = fast.cor(data[,x], data[,y], ndata)
       statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(ndata -3)
       p.value = pnorm(abs(statistic),
                   lower.tail = FALSE) * 2

    }#THEN
    # Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g") {

      statistic = mig.test(x, y, data, ndata, gsquare = TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Infomation (monte carlo permutation distribution)
    else if (test == "mc-mi") {

      perm.test = mc.test(data[,x], data[,y], ndata, samples = B, test = 1L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Pearson's X^2 test (monte carlo permutation distribution)
    else if (test == "mc-x2") {

      perm.test = mc.test(data[,x], data[,y], ndata, samples = B, test = 2L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Mutual Information for Gaussian Data (monte carlo permutation distribution)
    else if (test == "mc-mi-g") {

      statistic = mig.test(x, y, data, ndata, gsquare = TRUE)
      p.value = gmc.test(data[, x] , data[, y], B, test = 3L)

    }#THEN
    # Canonical (Linear) Correlation (monte carlo permutation distribution)
    else if (test == "mc-cor") {

      statistic = fast.cor(data[,x], data[,y], ndata)
      p.value = gmc.test(data[, x] , data[, y], B, test = 4L)

    }#THEN
    # Fisher's Z (monte carlo permutation distribution)
    else if (test == "mc-zf") {

      statistic = fast.cor(data[,x], data[,y], ndata)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(ndata -3)
      p.value = gmc.test(data[, x] , data[, y], B, test = 5L)

    }#THEN

  }#THEN
  else {

    # build the contingency table for discrete data only.
    if (test %in% available.discrete.tests) {

      # if there is only one parent, get it easy.
      if (length(sx) == 1)
        config = data[, sx]
      else
        config = configurations(data[, sx])

    }#THEN

    # Conditional Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi") {

      statistic = cmi.test(data[,x], data[,y], config, ndata, gsquare = TRUE)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Conditional Akaike Information Criterion-like test (binary, no p-value!)
    else if (test == "aict") {

      statistic = cmi.test(data[,x], data[,y], config, ndata, gsquare = FALSE) <
               (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config) / ndata
      p.value = as.integer(statistic)

    }#THEN
    # Shrinked Conditional Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi-sh") {

      statistic = shcmi.test(data[,x], data[,y], config, ndata, gsquare = TRUE)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (chi-square asymptotic distribution)
    else if (test == "x2") {

      statistic = cx2.test(data[,x], data[,y], config, ndata)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Canonical Partial Correlation (Student's t distribution)
    else if (test == "cor") {

      # define the degrees of freedom and check them.
      df = ndata - 2 - length(sx)
      if (df < 1)
        stop("trying to do a conditional independence test with zero degrees of freedom.")

      statistic = fast.pcor(x, y, sx, data, ndata)
      p.value = pt(abs(statistic * sqrt(df) / sqrt(1 - statistic^2)), df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      # define the degrees of freedom and check them.
      df = ndata - 3 - length(sx)
      if (df < 1)
        stop("trying to do a conditional independence test with zero degrees of freedom.")

      statistic = fast.pcor(x, y, sx, data, ndata)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(df)
      p.value = pnorm(abs(statistic), lower.tail = FALSE) * 2

    }#THEN
    # Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g") {

      statistic = cmig.test(x, y, sx, data, ndata, gsquare =TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Infomation (monte carlo permutation distribution)
    else if (test == "mc-mi") {

      perm.test = cmc.test(data[,x], data[,y], config, ndata, samples = B, test = 1L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Pearson's X^2 test (monte carlo permutation distribution)
    else if (test == "mc-x2") {

      perm.test = cmc.test(data[,x], data[,y], config, ndata, samples = B, test = 2L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Mutual Information for Gaussian Data (monte carlo permutation distribution)
    else if (test == "mc-mi-g") {

      statistic = cmig.test(x, y, sx, data, ndata, gsquare = TRUE)
      p.value = cgmc.test(x, y, sx, data, ndata, B, test = 3L)

    }#THEN
    # Canonical Partial Correlation (monte carlo permutation distribution)
    else if (test == "mc-cor") {

      statistic = fast.pcor(x, y, sx, data, ndata)
      p.value = cgmc.test(x, y, sx, data, ndata, B, test = 4L)

    }#THEN
    # Fisher's Z (monte carlo permutation distribution)
    else if (test == "mc-zf") {

      df = ndata - 3 - length(sx)
      statistic = fast.pcor(x, y, sx, data, ndata)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(df)
      p.value = cgmc.test(x, y, sx, data, ndata, B, test = 5L)

    }#THEN

  }#ELSE

  if (learning) {

    return(p.value)

  }#THEN
  else {

    if (test %in% c("aict")) {

      statistic = as.logical(statistic)
      p.value = NA

    }#THEN

    # build a valid object of class htest.
    result = structure(
        list(statistic = structure(statistic, names = test),
             p.value = p.value,
             method = test.labels[test],
             null.value = c(value = 0),
             alternative = ifelse(test %in% c("cor", "zf", "mc-cor", "mc-zf"),
               "two.sided", "greater"),
             data.name = paste(x, "~", y,
               ifelse(length(sx) > 0, "|", ""),
               paste(sx, collapse = " + "))
        ), class = "htest")

    if (!is.null(df))
      result$parameter = structure(df, names = "df")

    if (!is.null(B))
      result$parameter = structure(B, names = "Monte Carlo samples")

    return(result)

  }#ELSE

}#CONDITIONAL.TEST

# Mutual Information (discrete data)
mi.test = function(x, y, ndata, gsquare = TRUE) {

  s = .Call("mi",
      x = x,
      y = y,
      lx = nlevels(x),
      ly = nlevels(y),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#MI.TEST

# Shrinked Mutual Information (discrete data)
shmi.test = function(x, y, ndata, gsquare = TRUE) {

  s = .Call("shmi",
      x = x,
      y = y,
      lx = nlevels(x),
      ly = nlevels(y),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#SHMI.TEST

# Monte Carlo (discrete data)
mc.test = function(x, y, ndata, samples, test) {

  .Call("mcarlo",
        x = x,
        y = y,
        lx = nlevels(x),
        ly = nlevels(y),
        length = ndata,
        samples = samples,
        test = test,
        PACKAGE = "bnlearn")

}#MI.MC.TEST

# Conditional Mutual Information (discrete data)
cmi.test = function(x, y, z, ndata, gsquare = TRUE) {

  s = .Call("cmi",
      x = x,
      y = y,
      z = z,
      lx = nlevels(x),
      ly = nlevels(y),
      lz = nlevels(z),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#CMI.TEST

# Conditional Mutual Information (discrete data)
shcmi.test = function(x, y, z, ndata, gsquare = TRUE) {

  s = .Call("shcmi",
      x = x,
      y = y,
      z = z,
      lx = nlevels(x),
      ly = nlevels(y),
      lz = nlevels(z),
      length = ndata,
      PACKAGE = "bnlearn")

  ifelse(gsquare, 2 * ndata * s, s)

}#SHCMI.TEST

# Conditional Monte Carlo (discrete data)
cmc.test = function(x, y, z, ndata, samples, test) {

  .Call("cmcarlo",
        x = x,
        y = y,
        z = z,
        lx = nlevels(x),
        ly = nlevels(y),
        lz = nlevels(z),
        length = ndata,
        samples = samples,
        test = test,
        PACKAGE = "bnlearn")

}#CMI.MC.TEST

mig.test = function(x, y, data, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - fast.cor(data[, x], data[, y], ndata)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#MIG.TEST

# Conditional Mutual Information (gaussian data)
cmig.test = function(x, y, z, data, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - fast.pcor(x, y, z, data, ndata)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#CMIG.TEST

# Monte Carlo (gaussian data)
gmc.test = function(x, y, samples, test) {

  .Call("gauss_mcarlo",
        x = x,
        y = y,
        samples = samples,
        test = test,
        PACKAGE = "bnlearn")

}#GMC.TEST

# Conditional Monte Carlo (gaussian data)
cgmc.test = function(x, y, sx, data, ndata, samples, test) {

  .Call("gauss_cmcarlo",
        data = data[, c(x, y, sx), drop = FALSE],
        length = ndata,
        samples = samples,
        test = test,
        PACKAGE = "bnlearn")

}#CGMC.TEST

# Pearson's X^2 test (discrete data)
x2.test = function(x, y, ndata) {

  .Call("x2",
      x = x,
      y = y,
      lx = nlevels(x),
      ly = nlevels(y),
      length = ndata,
      PACKAGE = "bnlearn")

}#X2.TEST

# Pearson's Conditional X^2 test (discrete data)
cx2.test = function(x, y, z, ndata) {

  .Call("cx2",
      x = x,
      y = y,
      z = z,
      lx = nlevels(x),
      ly = nlevels(y),
      lz = nlevels(z),
      length = ndata,
      PACKAGE = "bnlearn")

}#CX2.TEST

# Fast implementation of the linear correlation coefficient.
fast.cor = function(x, y, ndata) {

  .Call("fast_cor",
        x = x,
        y = y,
        length = ndata,
        PACKAGE = "bnlearn")

}#FAST.COR

# Fast implementation of the partial correlation coefficient.
fast.pcor = function(x, y, sx, data, ndata) {

  .Call("fast_pcor",
        data = data[, c(x, y, sx), drop = FALSE],
        length = ndata,
        PACKAGE = "bnlearn")

}#FAST.PCOR

