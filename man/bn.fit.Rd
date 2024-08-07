\name{bn.fit}
\alias{bn.fit}
\alias{custom.fit}
\alias{bn.net}
\alias{$<-.bn.fit}
\title{Fit the parameters of a Bayesian network}
\description{

  Fit, assign or replace the parameters of a Bayesian network conditional on
  its structure.

}
\usage{
bn.fit(x, data, cluster, method, \dots, keep.fitted = TRUE,
  debug = FALSE)
custom.fit(x, dist, ordinal, debug = FALSE)
bn.net(x)
}
\arguments{
  \item{x}{an object of class \code{bn} (for \code{bn.fit()} and
    \code{custom.fit()}) or an object of class \code{bn.fit} (for
    \code{bn.net}).}
  \item{data}{a data frame containing the variables in the model.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{dist}{a named list, with element for each node of \code{x}. See below.}
  \item{method}{a character string, see below for details.}
  \item{...}{additional arguments for the parameter estimation procedure, see
    below.}
  \item{ordinal}{a vector of character strings, the labels of the discrete
    nodes which should be saved as ordinal random variables
    (\code{bn.fit.onode}) instead of unordered  factors (\code{bn.fit.dnode}).}
  \item{keep.fitted}{a boolean value. If \code{TRUE}, the object returned by
    \code{bn.fit} will contain fitted values and residuals for all Gaussian and
    conditional Gaussian nodes, and the configurations of the discrete parents
    for conditional Gaussian nodes.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{bn.fit()} fits the parameters of a Bayesian network given its structure
  and a data set; \code{bn.net} returns the structure underlying a fitted
  Bayesian network.

  \code{bn.fit()} accepts data with missing values encoded as \code{NA}. If the
  parameter estimation method was not specifically designed to deal with
  incomplete data, \code{bn.fit()} uses locally complete observations to fit the
  parameters of each local distribution.

  Available methods for \emph{discrete Bayesian networks} are:

  \itemize{

    \item \code{mle}: the maximum likelihood estimator for conditional
      probabilities.
    \item \code{bayes}: the classic Bayesian posterior estimator with a uniform
      prior matching that in the Bayesian Dirichlet equivalent (\code{bde})
      score.
    \item \code{hdir}: the hierarchical Dirichlet posterior estimator for
      related data sets from Azzimonti, Corani and Zaffalon (2019).
    \item \code{hard-em}: the Expectation-Maximization implementation of the
      estimators above.

  }

  Available methods for \emph{hybrid Bayesian networks} are:

  \itemize{

    \item \code{mle-g}: the maximum likelihood estimator for least squares
      regression models.
    \item \code{hard-em-g}: the Expectation-Maximization implementation of the
      estimators above.

  }

  Available methods for \emph{discrete Bayesian networks} are:

  \itemize{

    \item \code{mle-cg}: a combination of the maximum likelihood estimators
      \code{mle} and \code{mle-g}.
    \item \code{hard-em-cg}: the Expectation-Maximization implementation of the
      estimators above.

  }

  Additional arguments for the \code{bn.fit()} function:

  \itemize{

    \item \code{iss}: a numeric value, the imaginary sample size used by the
      \code{bayes} method to estimate the conditional probability tables
      associated with discrete nodes (see \code{\link{score}} for details).

    \item \code{replace.unidentifiable}: a boolean value. If \code{TRUE} and
      \code{method} one of \code{mle}, \code{mle-g} or \code{mle-cg},
      unidentifiable parameters are replaced by zeroes (in the case of
      regression coefficients and standard errors in Gaussian and conditional
      Gaussian nodes) or by uniform conditional probabilities (in discrete
      nodes).

      If \code{FALSE} (the default), the conditional probabilities in the local
      distributions of discrete nodes have a mximum likelihood estimate of
      \code{NaN} for all parents configurations that are not observed in
      \code{data}. Similarly, regression coefficients are set to \code{NA}
      if the linear regressions correspoding to the local distributions of
      continuous nodes are singular. Such missing values propagate to the
      results of functions such as \code{predict()}.

    \item \code{alpha0}: a positive number, the amount of information pooling
      between the related data sets in the \code{hdir} estimator.

    \item \code{group}: a character string, the label of the node with the
      grouping of the observations into the related data sets in the \code{hdir}
      estimator.

    \item \code{impute} and \code{impute.args}: a character string, the label of
      the imputation method (and its arguments) used by \code{hard-em},
      \code{hard-em-g} and \code{hard-em-cg} to complete the data in the
      \emph{expectation} step. The default method is the same as for
      \code{impute()}.

    \item \code{fit} and \code{fit.args}: a character string, the label of the
      parameter estimation method used by \code{hard-em}, \code{hard-em-g} and
      \code{hard-em-cg} to estimate the parameters in the \emph{maximization}
      step. The default method is the same as for \code{bn.fit()}.

    \item \code{loglik.threshold}: a non-negative numeric value, the minimum
      improvement threshold to continue iterating in \code{hard-em},
      \code{hard-em-g} and \code{hard-em-cg}. The threshold is defined as the
      relative likelihood improvement divided by the sample size of \code{data},
      and defaults to \code{1e-3}. Setting it to zero means that iterations only
      stop in case of negative improvement, which can happen due to stochastic
      noise if the imputation of the missing data uses approximate inference.

    \item \code{params.threshold}: a non-negative numeric value, the minimum
      maximum elative change in the parameter values to continue iterating in
      \code{hard-em}, \code{hard-em-g} and \code{hard-em-cg}. The threshold is
      defined as the maximum of the differences between parameter values divided
      scaled by the parameter value in the lastest iteration. The default value
      is \code{1e-3}.

    \item \code{max.iter}: a positive integer value, the maximum number of
      iterations in \code{hard-em}, \code{hard-em-g} and \code{hard-em-cg}. The
      default value is \code{5}.

    \item \code{start}: a \code{bn.fit} object, the fitted network used to
      initialize the \code{hard-em}, \code{hard-em-g} and \code{hard-em-cg}
      estimators. The default is to use the \code{bn.fit} object obtained from
      \code{x} with the default parameter estimator for the data, which will use
      locally complete data to fit the local distributions.

    \item \code{newdata}: a data frame, a separate set of data used to assess
      the convergence of the \code{hard-em}, \code{hard-em-g} and
      \code{hard-em-cg} estimators. The data in \code{data} are used by default
      for this purpose.

  }

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

  \code{custom.fit()} takes a set of user-specified distributions and their
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
  are no undirected arcs). See \code{\link{set.arc}} and \code{\link{cextend}}
  for two ways of manually setting the direction of one or more arcs.

  In the case of maximum likelihood estimators, \code{bn.fit()} produces
  \code{NA} parameter estimates for discrete and conditional Gaussian nodes when
  there are (discrete) parents configurations that are not observed in
  \code{data}. To avoid this either set \code{replace.unidentifiable} to
  \code{TRUE} or, in the case of discrete networks, use \code{method = "bayes"}.

}
\value{

  \code{bn.fit()} and \code{custom.fit()}returns an object of class
  \code{bn.fit}, \code{bn.net()} an object of class \code{bn}. See
  \code{\link{bn class}} and \code{\link{bn.fit class}} for details.

}
\seealso{
  \code{\link{bn.fit utilities}}, \code{\link{bn.fit plots}}.
}
\references{

  Azzimonti L, Corani G, Zaffalon M (2019). "Hierarchical Estimation of
    Parameters in Bayesian Networks".  \emph{Computational Statistics & Data
    Analysis}, \strong{137}:67--91.

}
\examples{
data(learning.test)

# learn the network structure.
cpdag = pc.stable(learning.test)
# set the direction of the only undirected arc, A - B.
dag = set.arc(cpdag, "A", "B")
# estimate the parameters of the Bayesian network.
fitted = bn.fit(dag, learning.test)
# replace the parameters of the node B.
new.cpt = matrix(c(0.1, 0.2, 0.3, 0.2, 0.5, 0.6, 0.7, 0.3, 0.1),
            byrow = TRUE, ncol = 3,
            dimnames = list(B = c("a", "b", "c"), A = c("a", "b", "c")))
fitted$B = as.table(new.cpt)
# the network structure is still the same.
all.equal(dag, bn.net(fitted))

# learn the network structure.
dag = hc(gaussian.test)
# estimate the parameters of the Bayesian network.
fitted = bn.fit(dag, gaussian.test)
# replace the parameters of the node F.
fitted$F = list(coef = c(1, 2, 3, 4, 5), sd = 3)
# set again the original parameters
fitted$F = lm(F ~ A + D + E + G, data = gaussian.test)

# discrete Bayesian network from expert knowledge.
dag = model2network("[A][B][C|A:B]")
cptA = matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
cptB = matrix(c(0.8, 0.2), ncol = 2, dimnames = list(NULL, c("GOOD", "BAD")))
cptC = c(0.5, 0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8)
dim(cptC) = c(2, 2, 2)
dimnames(cptC) = list("C" = c("TRUE", "FALSE"), "A" =  c("LOW", "HIGH"),
                   "B" = c("GOOD", "BAD"))
cfit = custom.fit(dag, dist = list(A = cptA, B = cptB, C = cptC))
# for ordinal nodes it is nearly the same.
cfit = custom.fit(dag, dist = list(A = cptA, B = cptB, C = cptC),
         ordinal = c("A", "B"))

# Gaussian Bayesian network from expert knowledge.
distA = list(coef = c("(Intercept)" = 2), sd = 1)
distB = list(coef = c("(Intercept)" = 1), sd = 1.5)
distC = list(coef = c("(Intercept)" = 0.5, "A" = 0.75, "B" = 1.32), sd = 0.4)
cfit = custom.fit(dag, dist = list(A = distA, B = distB, C = distC))

# conditional Gaussian Bayesian network from expert knowledge.
cptA = matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
distB = list(coef = c("(Intercept)" = 1), sd = 1.5)
distC = list(coef = matrix(c(1.2, 2.3, 3.4, 4.5), ncol = 2,
               dimnames = list(c("(Intercept)", "B"), NULL)),
          sd = c(0.3, 0.6))
cgfit = custom.fit(dag, dist = list(A = cptA, B = distB, C = distC))
}
\author{Marco Scutari}
\keyword{parameter learning}
