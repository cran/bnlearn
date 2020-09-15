\name{bn.cv}
\alias{bn.cv}
\alias{plot.bn.kcv}
\alias{plot.bn.kcv.list}
\alias{loss}
\title{Cross-validation for Bayesian networks}
\description{

  Perform a k-fold or hold-out cross-validation for a learning algorithm or a
  fixed network structure.

}
\usage{
bn.cv(data, bn, loss = NULL, ..., algorithm.args = list(),
  loss.args = list(), fit = "mle", fit.args = list(), method = "k-fold",
  cluster = NULL, debug = FALSE)

\method{plot}{bn.kcv}(x, ..., main, xlab, ylab, connect = FALSE)
\method{plot}{bn.kcv.list}(x, ..., main, xlab, ylab, connect = FALSE)

loss(x)
}
\arguments{
  \item{data}{a data frame containing the variables in the model.}
  \item{bn}{either a character string (the label of the learning algorithm to
    be applied to the training data in each iteration) or an object of class
    \code{bn} (a fixed network structure).}
  \item{loss}{a character string, the label of a loss function. If none is
    specified, the default loss function is the \emph{Classification Error}
    for Bayesian networks classifiers; otherwise, the \emph{Log-Likelihood
    Loss} for both discrete and continuous data sets. See below for
    additional details.}
  \item{algorithm.args}{a list of extra arguments to be passed to the learning
    algorithm.}
  \item{loss.args}{a list of extra arguments to be passed to the loss function
    specified by \code{loss}.}
  \item{fit}{a character string, the label of the method used to fit the
    parameters of the network. See \code{\link{bn.fit}} for details.}
  \item{fit.args}{additional arguments for the parameter estimation procedure,
    see again \code{\link{bn.fit}} for details.}
  \item{method}{a character string, either \code{k-fold}, \code{custom-folds}
    or \code{hold-out}. See below for details.}
  \item{cluster}{an optional cluster object from package \pkg{parallel}.}
  \item{debug}{a boolean value. If \code{TRUE} a lot of debugging output is
    printed; otherwise the function is completely silent.}
  \item{x}{an object of class \code{bn.kcv} or \code{bn.kcv.list} returned by
    \code{bn.cv()}.}
  \item{...}{additional objects of class \code{bn.kcv} or \code{bn.kcv.list}
    to plot alongside the first.}
  \item{main, xlab, ylab}{the title of the plot, an array of labels for the
    boxplot, the label for the y axis.}
  \item{connect}{a logical value. If \code{TRUE}, the medians points in the
    boxplots will be connected by a segmented line.}
}
\section{Cross-Validation Strategies}{

  The following cross-validation methods are implemented:
  \itemize{

    \item \emph{k-fold}: the \code{data} are split in \code{k} subsets of equal
      size. For each subset in turn, \code{bn} is fitted (and possibly learned
      as well) on the other \code{k - 1} subsets and the loss function is then
      computed using that subset. Loss estimates for each of the \code{k}
      subsets are then combined to give an overall loss for \code{data}.
    \item \emph{custom-folds}: the data are manually partitioned by the user
      into subsets, which are then used as in k-fold cross-validation. Subsets
      are not constrained to have the same size, and every observation must be
      assigned to one subset.
    \item \emph{hold-out}: \code{k} subsamples of size \code{m} are sampled
      independently without replacement from the \code{data}. For each subsample,
      \code{bn} is fitted (and possibly learned) on the remaining
      \code{m - nrow(data)} samples and the loss function is computed on the
      \code{m} observations in the subsample. The overall loss estimate is the
      average of the \code{k} loss estimates from the subsamples.
  }

  If cross-validation is used with multiple \code{runs}, the overall loss is the
  averge of the loss estimates from the different runs.

  To clarify, cross-validation methods accept the following optional arguments:
  \itemize{

    \item \code{k}: a positive integer number, the number of groups into which the
      data will be split (in k-fold cross-validation) or the number of times
      the data will be split in training and test samples (in hold-out
      cross-validation).
    \item \code{m}: a positive integer number, the size of the test set in
      hold-out cross-validation.
    \item \code{runs}: a positive integer number, the number of times
      k-fold or hold-out cross-validation will be run.
    \item \code{folds}: a list in which element corresponds to one fold and
      contains the indices for the observations that are included to that fold;
      or a list with an element for each run, in which each element is itself a
      list of the folds to be used for that run.

  }

}
\section{Loss Functions}{

  The following loss functions are implemented:

  \itemize{

    \item \emph{Log-Likelihood Loss} (\code{logl}): also known as \emph{negative
      entropy} or \emph{negentropy}, it is the negated expected log-likelihood
      of the test set for the Bayesian network fitted from the training set.
    \item \emph{Gaussian Log-Likelihood Loss} (\code{logl-g}): the negated
      expected log-likelihood for Gaussian Bayesian networks.
    \item \emph{Classification Error} (\code{pred}): the \emph{prediction error}
      for a single node in a discrete network. Frequentist predictions are used,
      so the values of the target node are predicted using only the information
      present in its local distribution (from its parents).
    \item \emph{Posterior Classification Error} (\code{pred-lw} and
      \code{pred-lw-cg}): similar to the above, but predictions are computed
      from an arbitrary set of nodes using likelihood weighting to obtain
      Bayesian posterior estimates. \code{pred-lw} applies to discrete Bayesian
      networks, \code{pred-lw-cg} to (discrete nodes in) hybrid networks.
    \item \emph{Exact Classification Error} (\code{pred-exact}): closed-form
      exact posterior predictions are available for Bayesian network classifiers.
    \item \emph{Predictive Correlation} (\code{cor}): the \emph{correlation}
      between the observed and the predicted values for a single node in a
      Gaussian Bayesian network.
    \item \emph{Posterior Predictive Correlation} (\code{cor-lw} and
      \code{cor-lw-cg}): similar to the above, but predictions are computed from
      an arbitrary set of nodes using likelihood weighting to obtain Bayesian
      posterior estimates. \code{cor-lw} applies to Gaussian networks and
      \code{cor-lw-cg} to (continuous nodes in) hybrid networks.
    \item \emph{Mean Squared Error} (\code{mse}): the \emph{mean squared error}
      between the observed and the predicted values for a single node in a
      Gaussian Bayesian network.
    \item \emph{Posterior Mean Squared Error} (\code{mse-lw} and
      \code{mse-lw-cg}): similar to the above, but predictions are computed from
      an arbitrary set of nodes using likelihood weighting to obtain Bayesian
      posterior estimates. \code{mse-lw} applies to Gaussian networks and
      \code{mse-lw-cg} to (continuous nodes in) hybrid networks.


  }

  Optional arguments that can be specified in \code{loss.args} are:

  \itemize{

    \item \code{target}: a character string, the label of target node for
      prediction in all loss functions but \code{logl}, \code{logl-g} and
      \code{logl-cg}.
    \item \code{from}: a vector of character strings, the labels of the nodes
      used to predict the \code{target} node in \code{pred-lw}, \code{pred-lw-cg},
      \code{cor-lw}, \code{cor-lw-cg}, \code{mse-lw} and \code{mse-lw-cg}. The
      default is to use all the other nodes in the network. Loss functions
      \code{pred}, \code{cor} and \code{mse} implicitly predict only from the
      parents of the \code{target} node.
    \item \code{n}: a positive integer, the number of particles used by
      likelihood weighting for \code{pred-lw}, \code{pred-lw-cg}, \code{cor-lw},
      \code{cor-lw-cg}, \code{mse-lw} and \code{mse-lw-cg}.
      The default value is \code{500}.

  }

  Note that if \code{bn} is a Bayesian network classifier, \code{pred} and
  \code{pred-lw} both give exact posterior predictions computed using the
  closed-form formulas for naive Bayes and TAN.

}
\section{Plotting Results from Cross-Validation}{

  Both plot methods accept any combination of objects of class \code{bn.kcv} or
  \code{bn.kcv.list} (the first as the \code{x} argument, the remaining as the
  \code{...} argument) and plot the respected expected loss values side by side.
  For a \code{bn.kcv} object, this mean a single point; for a \code{bn.kcv.list}
  object this means a boxplot.

}
\value{

  \code{bn.cv()} returns an object of class \code{bn.kcv.list} if \code{runs}
  is at least 2, an object of class \code{bn.kcv} if \code{runs} is equal to 1.

  \code{loss()} returns a numeric vector with a length equal to \code{runs}.

}
\references{

  Koller D, Friedman N (2009). \emph{Probabilistic Graphical Models: Principles
    and Techniques}. MIT Press.

}
\seealso{\code{\link{bn.boot}}, \code{\link{rbn}}, \code{\link{bn.kcv-class}}.}
\examples{
bn.cv(learning.test, 'hc', loss = "pred", loss.args = list(target = "F"))

folds = list(1:2000, 2001:3000, 3001:5000)
bn.cv(learning.test, 'hc', loss = "logl", method = "custom-folds",
  folds = folds)

xval = bn.cv(gaussian.test, 'mmhc', method = "hold-out",
         k = 5, m = 50, runs = 2)
xval
loss(xval)

\dontrun{
# comparing algorithms with multiple runs of cross-validation.
gaussian.subset = gaussian.test[1:50, ]
cv.gs = bn.cv(gaussian.subset, 'gs', runs = 10)
cv.iamb = bn.cv(gaussian.subset, 'iamb', runs = 10)
cv.inter = bn.cv(gaussian.subset, 'inter.iamb', runs = 10)
plot(cv.gs, cv.iamb, cv.inter,
  xlab = c("Grow-Shrink", "IAMB", "Inter-IAMB"), connect = TRUE)

# use custom folds.
folds = split(sample(nrow(gaussian.subset)), seq(5))
bn.cv(gaussian.subset, "hc", method = "custom-folds", folds = folds)

# multiple runs, with custom folds.
folds = replicate(5, split(sample(nrow(gaussian.subset)), seq(5)),
          simplify = FALSE)
bn.cv(gaussian.subset, "hc", method = "custom-folds", folds = folds)
}}
\author{Marco Scutari}
\keyword{resampling}
