
conditional.test = function(x, y, sx, data, test, learning = TRUE) {

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
    # Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1)

      if (obs.per.cell(x, y, data = data) >= 5) {

        statistic = mi.test(data[,x], data[,y], ndata, gsquare = TRUE)
        p.value = pchisq(statistic, df, lower.tail = FALSE)

      }#THEN
      else {

        statistic = 0
        p.value = 1

      }#ELSE

    }#THEN
    # Pearson's X^2 test (chi-square asymptotic distribution)
    else if (test == "x2") {

      statistic = x2.test(data[,x], data[,y], ndata)
      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1)
      p.value = pchisq(statistic, df, lower.tail = FALSE)

    }#THEN
    # Canonical (Linear) Correlation (Student's t distribution)
    else if (test == "cor") {

      statistic = cor(data[,x], data[,y])
      df = ndata - 2
      p.value = pt(abs((statistic * sqrt(ndata - 2) / sqrt(1 - statistic^2))),
                  df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

       statistic = cor(data[,x], data[,y])
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
    # Conditional Mutual Infomation a la FastIAMB (chi-square asymptotic distribution)
    else if (test == "fmi") {

      df = (nlevels(data[,x]) - 1) * (nlevels(data[,y]) - 1) * nlevels(config)

      if (obs.per.cell(x, y, config, data = data) >= 5) {

        statistic = cmi.test(data[,x], data[,y], config, ndata, gsquare = TRUE)
        p.value = pchisq(statistic, df, lower.tail = FALSE)

      }#THEN
      else {

        statistic = 0
        p.value = 1

      }#ELSE

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
        stop("tryng to do a conditional independence test with zero degrees of freedom.")

      statistic = pcor(c(x, y, sx), data)
      p.value = pt(abs(statistic * sqrt(df) / sqrt(1 - statistic^2)), df, lower.tail = FALSE) * 2

    }#THEN
    # Fisher's Z (asymptotic normal distribution)
    else if (test == "zf") {

      # define the degrees of freedom and check them.
      df = ndata - 3 - length(sx)
      if (df < 1)
        stop("tryng to do a conditional independence test with zero degrees of freedom.")

      statistic = pcor(c(x, y, sx), data)
      statistic = log((1 + statistic)/(1 - statistic))/2 * sqrt(df)
      p.value = pnorm(abs(statistic), lower.tail = FALSE) * 2

    }#THEN
    # Mutual Information for Gaussian Data (chi-square asymptotic distribution)
    else if (test == "mi-g") {

      statistic = cmig.test(x, y, sx, data, ndata, gsquare =TRUE)
      df = 1
      p.value = pchisq(statistic, df, lower.tail = FALSE)

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
             alternative = ifelse(test %in% c("cor", "zf"),
               "two.sided", "greater"),
             data.name = paste(x, "~", y,
               ifelse(length(sx) > 0, "|", ""),
               paste(sx, collapse = " + "))
        ), class = "htest")

    if (!is.null(df))
      result$parameter = structure(df, names = "df")

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

mig.test = function(x, y, data, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - cor(data[,x], data[,y])^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#MIG.TEST

# Conditional Mutual Information (gaussian data)
cmig.test = function(x, y, z, data, ndata, gsquare = TRUE) {

  s = - 0.5 * log(1 - pcor(c(x, y, z), data)^2)

  ifelse(gsquare, 2 * ndata * s, s)

}#CMIG.TEST

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

# Partial Canonical (Linear) Correlation.
# Original code was an internal copy of the pcor() function from
# package ggm; later adapted to use pseudoinverse() instead of
# solve() for greater reliability. Original copyright notice
# and licence were:
# Copyright Giovanni M. Marchetti, 2006.
# Released under "GPL version 2 or newer".
pcor = function (u, S) {

  k = pseudoinverse(cor(S[, u]))
  - k[1, 2] / sqrt(k[1, 1] * k[2, 2])

}#PCOR

# Moore-Penrose pseudoinverse matrix.
# Original code was an internal copy of the pseudoinverse()
# function from package corpcor version 1.4.6; removed unused
# stuff. Original copyright notice and licence were:
# Copyright Korbinian Strimmer 2003-2004.
# Released under "GPL version 2 or newer".
pseudoinverse = function (m) {

  msvd = positive.svd(m)

  if (length(msvd$d) == 0)
    array(0, dim(m)[2:1])
  else
    msvd$v %*% (1/msvd$d * t(msvd$u))

}#PSEUDOINVERSE

# Singular value decomposition that retains only positive singular values.
# Original code was an internal copy of the positive.svd()
# function from package corpcor version 1.4.6; removed unused
# stuff. Original copyright notice and licence were:
# Copyright Korbinian Strimmer 2003-2006.
# Released under "GPL version 2 or newer".
positive.svd = function(m) {

  s = svd(m)
  tol = max(dim(m)) * max(s$d) * .Machine$double.eps
  positive = s$d > tol

  list(d = s$d[positive], u = s$u[, positive, drop = FALSE],
      v = s$v[, positive, drop = FALSE])

}#POSITIVE.SVD

