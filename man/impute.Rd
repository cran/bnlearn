\name{impute}
\alias{impute}
\alias{predict.bn.fit}
\title{Predict or impute missing data from a Bayesian network}
\description{

  Impute missing values in a data set or predict a variable from a Bayesian
  network.

}
\usage{
\method{predict}{bn.fit}(object, node, data, method = "parents", ..., prob = FALSE,
  debug = FALSE)

impute(object, data, method, ..., debug = FALSE)
}
\arguments{
  \item{object}{an object of class \code{bn.fit} for \code{impute}; or an
    object of class \code{bn} or \code{bn.fit} for \code{predict}.}
  \item{data}{a data frame containing the data to be imputed. Complete
    observations will be ignored.}
  \item{node}{a character string, the label of a node.}
  \item{method}{a character string, the method used to impute the missing
    values or predict new ones. The default value is \code{parents}.}
  \item{\dots}{additional arguments for the imputation method. See below.}
  \item{prob}{a boolean value. If \code{TRUE} and \code{object} is a discrete
    network, the probabilities used for prediction are attached to the
    predicted values as an attribute called \code{prob}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
}
\details{

  \code{predict()} returns the predicted values for \code{node} given the data
  specified by \code{data} and the fitted network. Depending on the value of
  \code{method}, the predicted values are computed as follows.
  \itemize{
    \item \code{parents}: the predicted values are computed by plugging in
      the new values for the parents of \code{node} in the local probability
      distribution of \code{node} extracted from \code{fitted}.
    \item \code{bayes-lw}: the predicted values are computed by averaging
      likelihood weighting simulations performed using all the available nodes
      as evidence (obviously, with the exception of the node whose values we
      are predicting). The number of random samples which are averaged for each
      new observation is controlled by the \code{n} optional argument; the
      default is \code{500}. If the variable being predicted is discrete, the
      predicted level is that with the highest conditional probability. If the
      variable is continuous, the predicted value is the expected value of the
      conditional distribution. The variables that are used to compute the
      predicted values can be specified with the \code{from} optional argument;
      the default is to use all the relevant variables from the data. Note that
      the predicted values will differ in each call to \code{predict()} since
      this method is based on a stochastic simulation.
  }

  \code{impute()} is based on \code{predict()}, and can impute missing values
  with the same \code{methods} (\code{parents} and \code{bayes-lw}). The
  latter can take an additional argument \code{n} with the number of random
  samples which are averaged for each observation. As in \code{predict()},
  imputed values will differ in each call to \code{impute()} when \code{method}
  is set to \code{bayes-lw}.

  If \code{object} contains \code{NA} parameter estimates (because of
  unobserved discrete parents configurations in the data the parameters were
  learned from), \code{predict} will predict \code{NA}s when those parents
  configurations appear in \code{data}. See \code{\link{bn.fit}} for details on
  how to make sure \code{bn.fit} objects contain no \code{NA} parameter
  estimates.

}
\note{

  Ties in prediction are broken using \emph{Bayesian tie breaking}, i.e.
  sampling at random from the tied values. Therefore, setting the random seed
  is required to get reproducible results.

  Classifiers have a separate \code{predict()} method, see \link{naive.bayes}.

}
\value{

  \code{predict()} returns a numeric vector (for Gaussian and conditional
  Gaussian nodes), a factor (for categorical nodes) or an ordered factor (for
  ordinal nodes). If \code{prob = TRUE} and the network is discrete, the
  probabilities used for prediction are attached to the predicted values as
  an attribute called \code{prob}.

  \code{impute()} returns a data frame with the same structure as \code{data}.

}
\examples{
# missing data imputation.
with.missing.data = gaussian.test
with.missing.data[sample(nrow(with.missing.data), 500), "F"] = NA
fitted = bn.fit(model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]"),
           gaussian.test)
imputed = impute(fitted, with.missing.data)

# predicting a variable in the test set.
training = bn.fit(model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]"),
           gaussian.test[1:2000, ])
test = gaussian.test[2001:nrow(gaussian.test), ]
predicted = predict(training, node = "F", data = test)

# obtain the conditional probabilities for the values of a single variable
# given a subset of the rest, they are computed to determine the predicted
# values.
fitted = bn.fit(model2network("[A][C][F][B|A][D|A:C][E|B:F]"), learning.test)
evidence = data.frame(A = factor("a", levels = levels(learning.test$A)),
                      F = factor("b", levels = levels(learning.test$F)))
predicted = predict(fitted, "C", evidence,
              method = "bayes-lw", prob = TRUE)
attr(predicted, "prob")
}
\author{Marco Scutari}
\keyword{inference}
\keyword{missing data}
\keyword{data preprocessing}
