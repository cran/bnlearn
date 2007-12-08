# get all the subsets of a given size, even if either the initial set
# or the subset are empty (i.e. of size zero).
# Slightly modified version of the combinations() function from
# the gtools package. This is the copyroght notice:
## From email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
## dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
## Alex Ahgarin <datamanagement@email.com>.  Original version was
## named "subsets" and was Written by Bill Venables.  
# It's released under "LGPL 2.1".

subsets = function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {

  # allow empty subsets (i.e. subsets of empty sets).
  if ((n == 0) || (r == 0)) return(matrix(c(""),1,1))

  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 0)
    stop("bad value of n")

  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 0)
    stop("bad value of r")

  if (!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")

  if ((r > n) & repeats.allowed == FALSE)
    stop("r > n and repeats.allowed=FALSE")

  if (set) {

    v <- unique(sort(v))
    if (length(v) < n)
      stop("too few different elements")

  }#THEN

  v0 <- vector(mode(v), 0)
  if (repeats.allowed) {

    sub <- function(n, r, v) {

        if (r == 0)
            v0
        else if (r == 1)
            matrix(v, n, 1)
        else if (n == 1)
            matrix(v, 1, r)
        else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 1, r, v[-1]))

    }#SUB

  }#THEN
  else {

    sub <- function(n, r, v) {

        if (r == 0)
            v0
        else if (r == 1)
            matrix(v, n, 1)
        else if (r == n)
            matrix(v, 1, n)
        else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n - 1, r, v[-1]))

    }#SUB

  }#ELSE

  sub(n, r, v[1:n])

}#SUBSETS

# return the array whose size is smaller.
smaller = function(a, b) {

  if (length(a) < length(b))
    a
  else
    b

}#SMALLER

# build an array containing the configurations of the variables
# present in the 'data' argument.
configurations = function(data) {

  # this condition is there to avoid the remote possibility of an
  # integer overflow, which could happen if there are more than
  # MAXINT configurations. If you hit the 'else' branch, greetings,
  # you should not be here!
  if (prod(sapply(data, nlevels)) < 65536) {

    factor(.Call("cfg", data = data, PACKAGE = "bnlearn"))

  }#THEN
  else {

    warning("eeek! you have more than 65535 configurations!")
    factor(apply(data, 1, function(x) { 

      .Internal(paste(list(as.character(x)), sep = "", collapse = ":")) 

    }))

  }#ELSE

}#CONFIGURATIONS

