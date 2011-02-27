\name{naive.bayes}
\alias{naive.bayes}
\alias{predict.bn.naive}
\title{ Discrete naive Bayes classifiers }
\description{

  Create, fit and perform predictions with naive Bayes classifiers.

}
\usage{
naive.bayes(training, explanatory, data)
\method{predict}{bn.naive}(object, data, prior, ..., debug = FALSE)
}
\arguments{
  \item{training}{a character string, the label of the training variable.}
  \item{explanatory}{a vector of character strings, the labels of the
    explanatory variables.}
  \item{object}{an object of class \code{bn.naive}, either fitted or not.}
  \item{data}{a data frame containing the variables in the model.}
  \item{prior}{a numeric vector, the prior distribution for the training
    variable.}
 \item{\dots}{extra arguments from the generic method (currently ignored).}
 \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output
    is printed; otherwise the function is completely silent.}
}
\details{

  The \code{naive.bayes} functions creates the star-shaped Bayesian network
  form of a naive Bayes classifier; the training variable (the one holding
  the group each observation belongs to) is at the center of the star, and
  it has an outgoing arc for each explanatory variable.

  If \code{data} is specified, \code{explanatory} will be ignored and the
  labels of the explanatory variables will be extracted from the data.

  \code{predict} performs a supervised classification of the observations
  by assigning them to the group with the maximum posterior probability.

}
\note{

  Since \pkg{bnlearn} does not support networks containing both continuous
  and discrete variables, all variables in \code{data} must be discrete.

}
\value{

  \code{naive.bayes} returns an object of class \code{c("bn.naive", "bn")},
  which behaves like a normal \code{bn} object unless passed to \code{predict}.

  \code{predict} returns a factor with the same levels as the \code{training}
  variable from \code{data}.

}
\references{

  Borgelt C, Kruse R, Steinbrecher M (2009). \emph{Graphical Models: 
      Representations for Learning, Reasoning and Data Mining}.
      Wiley, 2nd edition.

}
\examples{
data(learning.test)
bn = naive.bayes("A", LETTERS[2:6])
pred = predict(bn, learning.test)
table(pred, learning.test[, "A"])
#     
# pred    a    b    c
#    a 1286  310  178
#    b  192  977  203
#    c  190  383 1281
# 
fitted = bn.fit(bn, learning.test, method = "bayes")
pred = predict(fitted, learning.test)
table(pred, learning.test[, "A"])
#     
# pred    a    b    c
#    a 1286  310  178
#    b  192  977  203
#    c  190  383 1281
#
}
\author{ Marco Scutari }
\keyword{multivariate}
