
conditional.test = function(x, y, sx, data, test, B, alpha = 1, learning = TRUE) {

  # update the test counter when performing structure learning.
  if (learning)
    increment.test.counter()

  sx = sx[sx != ""]
  ndata = nrow(data)
  df = NULL

  if (length(sx) == 0) {

    # all unconditional tests need this subsetting, do it here.
    datax = minimal.data.frame.column(data, x)
    datay = minimal.data.frame.column(data, y)

    # Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi") {

      par.test = mi.test(datax, datay, ndata, gsquare = TRUE)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Shrinked Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi-sh") {

      par.test = shmi.test(datax, datay, gsquare = TRUE)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (chi-square asymptotic distribution)
    else if (test == "x2") {

      par.test = x2.test(datax, datay)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Jonckheere-Terpstra Test (asymptotic normal distribution)
    else if (test == "jt") {

      statistic = jt(datax, datay)
      p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE)

    }#THEN
    # Canonical (Linear) Correlation (Student's t distribution)
    else if (test == "cor") {

      statistic = fast.cor(datax, datay, ndata)
      df = ndata - 2
      p.value = pt(abs((statistic * sqrt(ndata - 2) / sqrt(1 - statistic^2))),
                  df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      statistic = fast.cor(datax, datay, ndata)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(ndata -3)
      p.value = pnorm(abs(statistic), lower.tail = FALSE) * 2

    }#THEN
    # Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g") {

      statistic = mig.test(datax, datay, ndata, gsquare = TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Shrinked Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g-sh") {

      statistic = shmig.test(datax, datay, ndata, gsquare = TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Infomation (monte carlo permutation distribution)
    else if ((test == "mc-mi") || (test == "smc-mi")) {

      perm.test = mc.test(datax, datay, samples = B,
                    alpha = ifelse(test == "smc-mi", alpha, 1), test = 1L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Mutual Infomation (semiparametric)
    else if (test == "sp-mi") {

      perm.test = mc.test(datax, datay, samples = B,
                    alpha = ifelse(test == "smc-mi", alpha, 1), test = 6L)
      statistic = perm.test[1]
      df = perm.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (monte carlo permutation distribution)
    else if ((test == "mc-x2") || (test == "smc-x2")) {

      perm.test = mc.test(datax, datay, samples = B,
                    alpha = ifelse(test == "smc-x2", alpha, 1), test = 2L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Pearson's X^2 test (semiparametric)
    else if (test == "sp-x2") {

      perm.test = mc.test(datax, datay, samples = B,
                    alpha = ifelse(test == "smc-mi", alpha, 1), test = 7L)
      statistic = perm.test[1]
      df = perm.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Information for Gaussian Data (monte carlo permutation distribution)
    else if ((test == "mc-mi-g") || (test == "smc-mi-g")) {

      statistic = mig.test(datax, datay, ndata, gsquare = TRUE)
      p.value = gmc.test(datax, datay, samples = B,
                  alpha = ifelse(test == "smc-mi-g", alpha, 1), test = 3L)

    }#THEN
    # Canonical (Linear) Correlation (monte carlo permutation distribution)
    else if ((test == "mc-cor") || (test == "smc-cor")) {

      statistic = fast.cor(datax, datay, ndata)
      p.value = gmc.test(datax, datay, samples = B,
                  alpha = ifelse(test == "smc-cor", alpha, 1), test = 4L)

    }#THEN
    # Fisher's Z (monte carlo permutation distribution)
    else if ((test == "mc-zf") || (test == "smc-zf")) {

      statistic = fast.cor(datax, datay, ndata)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(ndata -3)
      p.value = gmc.test(datax, datay, samples = B,
                  alpha = ifelse(test == "smc-zf", alpha, 1), test = 5L)

    }#THEN
    # Jonckheere-Terpstra test (monte carlo permutation distribution)
    else if ((test == "mc-jt") || (test == "smc-jt")) {

      perm.test = mc.test(datax, datay, samples = B,
                    alpha = ifelse(test == "smc-jt", alpha, 1), test = 8L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN

  }#THEN
  else {

    # build the contingency table only for discrete data.
    if (test %in% c(available.discrete.tests, available.ordinal.tests)) {

      # if there is only one parent, get it easy.
      if (length(sx) == 1)
        config = minimal.data.frame.column(data, sx)
      else
        config = configurations(minimal.data.frame.column(data, sx))

    }#THEN

    # Conditional Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      par.test = cmi.test(datax, datay, config, ndata, gsquare = TRUE)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Shrinked Conditional Mutual Infomation (chi-square asymptotic distribution)
    if (test == "mi-sh") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      par.test = shcmi.test(datax, datay, config, gsquare = TRUE)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (chi-square asymptotic distribution)
    else if (test == "x2") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      par.test = cx2.test(datax, datay, config)
      statistic = par.test[1]
      df = par.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Jonckheere-Terpstra Test (asymptotic normal distribution)
    else if (test == "jt") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      statistic = cjt(datax, datay, config)
      p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE)

    }#THEN
    # Canonical Partial Correlation (Student's t distribution)
    else if (test == "cor") {

      # define the degrees of freedom and check them.
      df = ndata - 2 - length(sx)
      if (df < 1)
        stop("trying to do a conditional independence test with zero degrees of freedom.")

      statistic = fast.pcor(x, y, sx, data, ndata, strict = !learning)
      p.value = pt(abs(statistic * sqrt(df) / sqrt(1 - statistic^2)), df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      # define the degrees of freedom and check them.
      df = ndata - 3 - length(sx)
      if (df < 1)
        stop("trying to do a conditional independence test with zero degrees of freedom.")

      statistic = fast.pcor(x, y, sx, data, ndata, strict = !learning)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(df)
      p.value = pnorm(abs(statistic), lower.tail = FALSE) * 2

    }#THEN
    # Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g") {

      statistic = cmig.test(x, y, sx, data, ndata, strict = !learning, gsquare = TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Shrinked Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g-sh") {

      statistic = shcmig.test(x, y, sx, data, ndata, gsquare = TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Infomation (monte carlo permutation distribution)
    else if ((test == "mc-mi") || (test == "smc-mi")) {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      perm.test = cmc.test(datax, datay, config, samples = B,
                    alpha = ifelse(test == "smc-mi", alpha, 1), test = 1L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Mutual Infomation (semiparametric)
    else if (test == "sp-mi") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      perm.test = cmc.test(datax, datay, config, samples = B,
                    alpha = ifelse(test == "smc-mi", alpha, 1), test = 6L)
      statistic = perm.test[1]
      df = perm.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Pearson's X^2 test (monte carlo permutation distribution)
    else if ((test == "mc-x2") || (test == "smc-x2")) {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      perm.test = cmc.test(datax, datay, config, samples = B,
                    alpha = ifelse(test == "smc-x2", alpha, 1), test = 2L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN
    # Pearson's X^2 test (semiparametric)
    else if (test == "sp-x2") {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      perm.test = cmc.test(datax, datay, config, samples = B,
                    alpha = ifelse(test == "smc-x2", alpha, 1), test = 7L)
      statistic = perm.test[1]
      df = perm.test[2]
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Mutual Information for Gaussian Data (monte carlo permutation distribution)
    else if ((test == "mc-mi-g") || (test == "smc-mi-g")) {

      statistic = cmig.test(x, y, sx, data, ndata, strict = !learning, gsquare = TRUE)
      p.value = cgmc.test(x, y, sx, data, ndata, samples = B,
                  alpha = ifelse(test == "smc-mi-g", alpha, 1), test = 3L)

    }#THEN
    # Canonical Partial Correlation (monte carlo permutation distribution)
    else if ((test == "mc-cor") || (test == "smc-cor")) {

      statistic = fast.pcor(x, y, sx, data, ndata, strict = !learning)
      p.value = cgmc.test(x, y, sx, data, ndata, samples = B,
                  alpha = ifelse(test == "smc-cor", alpha, 1), test = 4L)

    }#THEN
    # Fisher's Z (monte carlo permutation distribution)
    else if ((test == "mc-zf") || (test == "smc-zf")) {

      df = ndata - 3 - length(sx)
      statistic = fast.pcor(x, y, sx, data, ndata, strict = !learning)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(df)
      p.value = cgmc.test(x, y, sx, data, ndata, samples = B,
                  alpha = ifelse(test == "smc-zf", alpha, 1), test = 5L)

    }#THEN
    # Jonckheere-Terpstra test (monte carlo permutation distribution)
    else if ((test == "mc-jt") || (test == "smc-jt")) {

      datax = minimal.data.frame.column(data, x)
      datay = minimal.data.frame.column(data, y)

      perm.test = cmc.test(datax, datay, config, samples = B,
                    alpha = ifelse(test == "smc-jt", alpha, 1), test = 8L)
      statistic = perm.test[1]
      p.value = perm.test[2]

    }#THEN

  }#ELSE

  if (learning) {

    return(p.value)

  }#THEN
  else {

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

    # degrees of freedom.
    if (!is.null(df))
      result$parameter = structure(df, names = "df")
    # number of permutations.
    if (!is.null(B))
      result$parameter = c(result$parameter, structure(B, names = "Monte Carlo samples"))

    return(result)

  }#ELSE

}#CONDITIONAL.TEST

# Mutual Information (discrete data)
mi.test = function(x, y, ndata, gsquare = TRUE) {

  s = .Call("mi",
      x = x,
      y = y,
      gsquare = gsquare)

}#MI.TEST

# Shrinked Mutual Information (discrete data)
shmi.test = function(x, y, gsquare = TRUE) {

  s = .Call("shmi",
      x = x,
      y = y,
      gsquare = gsquare)

}#SHMI.TEST

# Monte Carlo (discrete data)
mc.test = function(x, y, samples, alpha, test) {

  .Call("mcarlo",
        x = x,
        y = y,
        samples = samples,
        test = test,
        alpha = alpha)

}#MC.TEST

# Conditional Mutual Information (discrete data)
cmi.test = function(x, y, z, ndata, gsquare = TRUE) {

  s = .Call("cmi",
      x = x,
      y = y,
      z = z,
      gsquare = gsquare)

}#CMI.TEST

# Shrinked Conditional Mutual Information (discrete data)
shcmi.test = function(x, y, z, gsquare = TRUE) {

  s = .Call("shcmi",
      x = x,
      y = y,
      z = z,
      gsquare = gsquare)

}#SHCMI.TEST

# Conditional Monte Carlo (discrete data)
cmc.test = function(x, y, z, samples, alpha, test) {

  .Call("cmcarlo",
        x = x,
        y = y,
        z = z,
        samples = samples,
        test = test,
        alpha = alpha)

}#CMC.TEST

# Mutual Information (gaussian data)
mig.test = function(x, y, ndata, gsquare = TRUE) {

  s = .Call("mig",
      x = x,
      y = y,
      length = ndata)

  ifelse(gsquare, 2 * ndata * s, s)

}#MIG.TEST

# Shrinked Mutual Information (gaussian data)
shmig.test = function(x, y, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - fast.shcor(x, y, ndata)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#SHMIG.TEST

# Conditional Mutual Information (gaussian data)
cmig.test = function(x, y, z, data, ndata, strict, gsquare = TRUE) {

  s = - 0.5 * log(1 - fast.pcor(x, y, z, data, ndata, strict)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#CMIG.TEST

# Shrinked Conditional Mutual Information (gaussian data)
shcmig.test = function(x, y, z, data, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - fast.shpcor(x, y, z, data, ndata, strict = TRUE)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#SHCMIG.TEST

# Monte Carlo (gaussian data)
gmc.test = function(x, y, samples, alpha, test) {

  .Call("gauss_mcarlo",
        x = x,
        y = y,
        samples = samples,
        test = test,
        alpha = alpha)

}#GMC.TEST

# Conditional Monte Carlo (gaussian data)
cgmc.test = function(x, y, sx, data, ndata, samples, alpha, test) {

  .Call("gauss_cmcarlo",
        data = minimal.data.frame.column(data, c(x, y, sx)),
        length = ndata,
        samples = samples,
        test = test,
        alpha = alpha)

}#CGMC.TEST

# Pearson's X^2 test (discrete data)
x2.test = function(x, y) {

  .Call("x2",
      x = x,
      y = y)

}#X2.TEST

# Pearson's Conditional X^2 test (discrete data)
cx2.test = function(x, y, z) {

  .Call("cx2",
      x = x,
      y = y,
      z = z)

}#CX2.TEST

# Fast implementation of the linear correlation coefficient.
fast.cor = function(x, y, ndata) {

  .Call("fast_cor",
        x = x,
        y = y,
        length = ndata)

}#FAST.COR

# Fast implementation of the partial correlation coefficient.
fast.pcor = function(x, y, sx, data, ndata, strict) {

  .Call("fast_pcor",
        data = minimal.data.frame.column(data, c(x, y, sx)),
        length = ndata,
        shrinkage = FALSE,
        strict = strict)

}#FAST.PCOR

# Fast implementation of the shrinked linear correlation coefficient.
fast.shcor = function(x, y, ndata) {

  .Call("fast_shcor",
        x = x,
        y = y,
        length = ndata)

}#FAST.SHCOR

# Fast implementation of the shrinked partial correlation coefficient.
fast.shpcor = function(x, y, sx, data, ndata, strict) {

  .Call("fast_pcor",
        data = minimal.data.frame.column(data, c(x, y, sx)),
        length = ndata,
        shrinkage = TRUE,
        strict = strict)

}#FAST.SHPCOR

# Jonckheere-Terpstra test.
jt = function(x, y) {

  .Call("jt",
        x = x,
        y = y)

}#JT

# Conditional Jonckheere-Terpstra test.
cjt = function(x, y, z) {

  .Call("cjt",
        x = x,
        y = y,
        z = z)

}#CJT
