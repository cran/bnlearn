\name{naive.bayes}
\alias{naive.bayes}
\alias{tree.bayes}
\alias{predict.bn.naive}
\alias{predict.bn.tan}
\title{ Naive Bayes classifiers }
\description{

  Create, fit and perform predictions with naive Bayes and 
  Tree-Augmented naive Bayes (TAN) classifiers.

}
\usage{
naive.bayes(training, explanatory, data)
\method{predict}{bn.naive}(object, data, prior, ..., debug = FALSE)

tree.bayes(x, training, explanatory, whitelist = NULL, blacklist = NULL,
  mi = NULL, root = NULL, debug = FALSE)
\method{predict}{bn.tan}(object, data, prior, ..., debug = FALSE)
}
\arguments{
  \item{training}{a character string, the label of the training variable.}
  \item{explanatory}{a vector of character strings, the labels of the
    explanatory variables.}
  \item{object}{an object of class \code{bn.naive}, either fitted or not.}
  \item{x, data}{a data frame containing the variables in the model, which
    must all be factors.}
  \item{prior}{a numeric vector, the prior distribution for the training
    variable. It is automatically normalized if not already so.}
  \item{whitelist}{a data frame with two columns (optionally labeled
    "from" and "to"), containing a set of arcs to be included in the
     graph. }
  \item{blacklist}{a data frame with two columns (optionally labeled
    "from" and "to"), containing a set of arcs not to be included in
     the graph. }
  \item{mi}{a character string, the estimator used for the mutual
     information coefficients for the Chow-Liu algorithm in TAN. Possible
     values are \code{mi} (discrete mutual information) and \code{mi-g}
     (Gaussian mutual information).}
 \item{root}{a character string, the label of the explanatory variable to
     bre used as the root of the tree in the TAN classifier.}
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

  Ties in prediction are broken using \emph{Bayesian tie breaking}, i.e.
  sampling at random from the tied values. Therefore, setting the random
  seed is required to get reproducible results.

}
\value{

  \code{naive.bayes} returns an object of class \code{c("bn.naive", "bn")},
  which behaves like a normal \code{bn} object unless passed to \code{predict}.
  \code{tree.bayes} returns an object of class \code{c("bn.tan", "bn")},
  which again behaves like a normal \code{bn} object unless passed to 
  \code{predict}.

  \code{predict} returns a factor with the same levels as the \code{training}
  variable from \code{data}.

}
\references{

  Borgelt C, Kruse R, Steinbrecher M (2009). \emph{Graphical Models:
      Representations for Learning, Reasoning and Data Mining}.
      Wiley, 2nd edition.

  Friedman N, Geiger D, Goldszmidt M (1997). "Bayesian Network Classifiers".
      \emph{Machine Learning}, \strong{29}(2--3), 131--163. 

}
\examples{
data(learning.test)
bn = naive.bayes("A", LETTERS[2:6])
pred = predict(bn, learning.test)
table(pred, learning.test[, "A"])
# pred    a    b    c
#    a 1286  310  178
#    b  192  977  203
#    c  190  383 1281

tan = tree.bayes(learning.test, "A")

fitted = bn.fit(tan, learning.test, method = "bayes")
pred = predict(fitted, learning.test)
table(pred, learning.test[, "A"])
# pred    a    b    c
#    a 1300  180  184
#    b  206 1342  239
#    c  162  148 1239
}
\author{ Marco Scutari }
\keyword{multivariate}
