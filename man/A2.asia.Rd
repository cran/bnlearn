\name{asia}
\docType{data}
\alias{asia}
\title{ Asia (synthetic) data set by Lauritzen and Spiegelhalter}
\description{

  Small synthetic data set from Lauritzen and Spiegelhalter (1988)
  about lung diseases (tuberculosis, lung cancer or bronchitis)
  and visits to Asia.

}
\usage{
  data(asia)
}
\format{

  The \code{asia} data set contains the following variables:
  \itemize{

    \item \code{D} (\emph{dyspnoea}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{T} (\emph{tuberculosis}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{L} (\emph{lung cancer}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{B} (\emph{bronchitis}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{A} (\emph{visit to Asia}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{S} (\emph{smoking}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{X} (\emph{chest X-ray}), a two-level factor with levels
            \code{yes} and \code{no}.
    \item \code{E} (\emph{tuberculosis versus lung cancer/bronchitis}), a
            two-level factor with levels \code{yes} and \code{no}.

  }

}
\note{

  Standard learning algorithms are not able to recover the true
  structure of the network because of the presence of a node (\code{E})
  with conditional probabilities equal to both 0 and 1.

}
\source{

  S. Lauritzen and D. Spiegelhalter (1988). Local computation with
    probabilities on graphical structures and their application to expert
    system. Journal of the Royal Statistics Society - B Series, 50(2),
    pages 157--192.

}
\examples{

## The modelstring() of this data set is:
# [A][S][T|A][L|S][B|S][D|B][E|T:L][X|E]

# these are the R commands used to generate this data set.
\dontrun{
a = sample(c("yes", "no"), 5000, prob = c(0.01, 0.99), replace = TRUE)
s = sample(c("yes", "no"), 5000, prob = c(0.50, 0.50), replace = TRUE)

t = a
t[t == "yes"] = sample(c("yes", "no"), length(which(t == "yes")),
                prob = c(0.05, 0.95), replace = TRUE)
t[t == "no"] = sample(c("yes", "no"), length(which(t == "no")),
                prob = c(0.01, 0.99), replace = TRUE)

l = s
l[l == "yes"] = sample(c("yes", "no"), length(which(l == "yes")),
                prob = c(0.10, 0.90), replace = TRUE)
l[l == "no"] = sample(c("yes", "no"), length(which(l == "no")),
                prob = c(0.01, 0.99), replace = TRUE)

b = s
b[b == "yes"] = sample(c("yes", "no"), length(which(b == "yes")),
                prob = c(0.60, 0.40), replace = TRUE)
b[b == "no"] = sample(c("yes", "no"), length(which(b == "no")),
                prob = c(0.30, 0.70), replace = TRUE)

e = apply(cbind(l,t), 1, paste, collapse= ":")
e[e == "yes:yes"] = "yes"
e[e == "yes:no"] = "yes"
e[e == "no:yes"] = "yes"
e[e == "no:no"] = "no"

x = e
x[x == "yes"] = sample(c("yes", "no"), length(which(x == "yes")),
                prob = c(0.98, 0.02), replace = TRUE)
x[x == "no"] = sample(c("yes", "no"), length(which(x == "no")),
                prob = c(0.05, 0.95), replace = TRUE)

d = apply(cbind(e,b), 1, paste, collapse= ":")
d[d == "yes:yes"] = sample(c("yes", "no"), length(which(d == "yes:yes")),
                    prob = c(0.90, 0.10), replace = TRUE)
d[d == "yes:no"] = sample(c("yes", "no"), length(which(d == "yes:no")),
                    prob = c(0.70, 0.30), replace = TRUE)
d[d == "no:yes"] = sample(c("yes", "no"), length(which(d == "no:yes")),
                   prob = c(0.80, 0.20), replace = TRUE)
d[d == "no:no"] = sample(c("yes", "no"), length(which(d == "no:no")),
                  prob = c(0.10, 0.90), replace = TRUE)

data.frame(A = a, S = s, T = t, L = l, B = b, E = e, X = x, D = d)
}

}
\keyword{datasets}
