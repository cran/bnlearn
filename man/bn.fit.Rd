\name{bn.fit}
\alias{bn.fit}
\alias{custom.fit}
\alias{bn.net}
\alias{$<-.bn.fit}
\title{Fit the parameters of a Bayesian network}
\description{

  Fit the parameters of a Bayesian network conditional on its structure.

}
\usage{
bn.fit(x, data, method = "mle", \dots, debug = FALSE)
custom.fit(x, dist, ordinal)
bn.net(x, debug = FALSE)
}
\arguments{
  \item{x}{an object of class \code{bn} (for \code{bn.fit} and \code{custom.fit})
    or an object of class \code{bn.fit} (for \code{bn.net}).}
  \item{data}{a data frame containing the variables in the model.}
  \item{dist}{a named list, with element for each node of \code{x}. See below.}
  \item{method}{a character string, either \code{mle} for \emph{Maximum
    Likelihood parameter estimation} or \code{bayes} for \emph{Bayesian
    parameter estimation} (currently implemented only for discrete data).}
  \item{...}{additional arguments for the parameter estimation prcoedure, see
    below.}
  \item{ordinal}{a vector of character strings, the labels of the discrete
    nodes which should be saved as ordinal random variables
    (\code{bn.fit.onode}) instead of unordered  factors (\code{bn.fit.dnode}).}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{bn.fit} fits the parameters of a Bayesian network given its structure
  and a data set; \code{bn.net} returns the structure underlying a fitted 
  Bayesian network.

  An in-place replacement method is available to change the parameters of each
  node in a \code{bn.fit} object; see the examples for discrete, continuous and
  hybrid networks below. For a discrete node (class \code{bn.fit.dnode} or
  \code{bn.fit.onode}), the new parameters must be in a \code{table} object.
  For a Gaussian node (class \code{bn.fit.gnode}), the new parameters can be
  defined either by an \code{lm}, \code{glm} or \code{pensim} object (the
  latter is from the \code{penalized} package) or in a list with elements named
  \code{coef}, \code{sd} and optionally \code{fitted} and \code{resid}. For
  a conditional Gaussian node (class \code{bn.fit.cgnode}), the new parameters
  can be defined by a list with elements named \code{coef}, \code{sd} and
  optionally \code{fitted}, \code{resid} and \code{configs}. In both cases
  \code{coef} should contain the new regression coefficients, \code{sd} the
  standard deviation of the residuals, \code{fitted} the fitted values and
  \code{resid} the residuals. \code{configs} should contain the configurations
  if the discrete parents of the conditional Gaussian node, stored as a factor.

  \code{custom.fit} takes a set of user-specified distributions and their
  parameters and uses them to build a \code{bn.fit} object. Its purpose is to
  specify a Bayesian network (complete with the parameters, not only the 
  structure) using knowledge from experts in the field instead of learning it
  from a data set. The distributions must be passed to the function in a list,
  with elements named after the nodes of the network structure \code{x}. Each
  element of the list must be in one of the formats described above for
  in-place replacement.

}
\note{

  Due to the way Bayesian networks are defined it is possible to estimate their
  parameters only if the network structure is completely directed (i.e. there
  are no undirected arcs). See \code{\link{set.arc}} and \code{\link{pdag2dag}}
  for two ways of manually setting the direction of one or more arcs.

  The only supported additional parameter is the imaginary sample size
  (\code{iss}) for the Dirichlet posterior distribution of discrete networks
  (see \code{\link{score}} for details).

  The conditional probabilities in the local distributions of discrete nodes
  have a mximum likelihood estimate of \code{NaN} for all parents configurations
  that are not observed in \code{data}. Such missing values propagate to the
  results of functions such as \code{predict}. Bayesian posterior estimates do
  not have this problem, and are very close to the corresponding maximum
  likelihood estimates when \code{iss} is small.

}
\value{

  \code{bn.fit} returns an object of class \code{bn.fit}, \code{bn.net}
  an object of class \code{bn}. See \code{\link{bn class}} and
  \code{\link{bn.fit class}} for details.

}
\seealso{
  \code{\link{bn.fit utilities}}, \code{\link{bn.fit plots}}.
}
\examples{
data(learning.test)

# learn the network structure.
res = gs(learning.test)
# set the direction of the only undirected arc, A - B.
res = set.arc(res, "A", "B")
# estimate the parameters of the Bayesian network.
fitted = bn.fit(res, learning.test)
# replace the parameters of the node B.
new.cpt = matrix(c(0.1, 0.2, 0.3, 0.2, 0.5, 0.6, 0.7, 0.3, 0.1),
            byrow = TRUE, ncol = 3,
            dimnames = list(B = c("a", "b", "c"), A = c("a", "b", "c")))
fitted$B = as.table(new.cpt)
# the network structure is still the same.
all.equal(res, bn.net(fitted))

# learn the network structure.
res = hc(gaussian.test)
# estimate the parameters of the Bayesian network.
fitted = bn.fit(res, gaussian.test)
# replace the parameters of the node F.
fitted$F = list(coef = c(1, 2, 3, 4, 5), sd = 3)
# set again the original parameters
fitted$F = lm(F ~ A + D + E + G, data = gaussian.test)

# discrete Bayesian network from expert knowledge.
net = model2network("[A][B][C|A:B]")
cptA = matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
cptB = matrix(c(0.8, 0.2), ncol = 2, dimnames = list(NULL, c("GOOD", "BAD")))
cptC = c(0.5, 0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8)
dim(cptC) = c(2, 2, 2)
dimnames(cptC) = list("C" = c("TRUE", "FALSE"), "A" =  c("LOW", "HIGH"),
                   "B" = c("GOOD", "BAD"))
cfit = custom.fit(net, dist = list(A = cptA, B = cptB, C = cptC))
# for ordinal nodes it is nearly the same.
cfit = custom.fit(net, dist = list(A = cptA, B = cptB, C = cptC),
         ordinal = c("A", "B"))

# Gaussian Bayesian network from expert knowledge.
distA = list(coef = c("(Intercept)" = 2), sd = 1)
distB = list(coef = c("(Intercept)" = 1), sd = 1.5)
distC = list(coef = c("(Intercept)" = 0.5, "A" = 0.75, "B" = 1.32), sd = 0.4)
cfit = custom.fit(net, dist = list(A = distA, B = distB, C = distC))

# conditional Gaussian Bayesian network from expert knowledge.
cptA = matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
distB = list(coef = c("(Intercept)" = 1), sd = 1.5)
distC = list(coef = matrix(c(1.2, 2.3, 3.4, 4.5), ncol = 2,
               dimnames = list(c("(Intercept)", "B"), NULL)),
          sd = c(0.3, 0.6))
cgfit = custom.fit(net, dist = list(A = cptA, B = distB, C = distC))
}
\author{ Marco Scutari }
\keyword{multivariate}
