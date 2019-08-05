\name{independence-tests}
\alias{independence-tests}
\alias{independence tests}
\title{Conditional independence tests}
\description{

  Overview of the conditional independence tests implemented in \pkg{bnlearn},
  with the respective reference publications.

}
\details{

  Unless otherwise noted, the reference publication for conditional independence
  tests is:

  Edwards DI (2000). \emph{Introduction to Graphical Modelling}. Springer, 2nd
    edition.

  Additionally for continuous permutation tests:

  Legendre P (2000). "Comparison of Permutation Methods for the Partial
    Correlation and Partial Mantel Tests". \emph{Journal of Statistical
    Computation and Simulation}, \strong{67}:37--73.

  and for semiparametric discrete tests:

  Tsamardinos I, Borboudakis G (2010). "Permutation Testing Improves Bayesian
    Network Learning". \emph{Machine Learning and Knowledge Discovery in
    Databases}, 322--337.

  Available conditional independence tests (and the respective labels) for
  \emph{discrete Bayesian networks} (categorical variables) are:

  \itemize{

    \item \emph{mutual information}: an information-theoretic distance measure.
	  It's proportional to the log-likelihood ratio (they differ by a
	  \eqn{2n}{2n} factor) and is related to the deviance of the tested models.
	  The asymptotic \eqn{\chi^2}{chi-square} test (\code{mi} and
	  \code{mi-adf}, with adjusted degrees of freedom), the Monte Carlo
	  permutation test (\code{mc-mi}), the sequential Monte Carlo permutation
	  test (\code{smc-mi}), and the semiparametric test (\code{sp-mi}) are
	  implemented.

    \item \emph{shrinkage estimator} for the \emph{mutual information}
      (\code{mi-sh}): an improved asymptotic \eqn{\chi^2}{chi-square} test
      based on the James-Stein estimator for the mutual information.

      Hausser J, Strimmer K (2009). "Entropy inference and the James-Stein
        estimator, with application to nonlinear gene association networks".
        \emph{Statistical Applications in Genetics and Molecular Biology},
        \strong{10}:1469--1484.

    \item \emph{Pearson's \eqn{X^2}{X^2}}: the classical Pearson's
	  \eqn{X^2}{X^2} test for contingency tables. The asymptotic
	  \eqn{\chi^2}{chi-square} test (\code{x2} and \code{x2-adf}, with adjusted
	  degrees of freedom), the Monte Carlo permutation test (\code{mc-x2}), the
	  sequential Monte Carlo permutation test (\code{smc-x2}) and
	  semiparametric test (\code{sp-x2}) are implemented.

  }

  Available conditional independence tests (and the respective labels) for
  \emph{discrete Bayesian networks} (ordered factors) are:

  \itemize{

    \item \emph{Jonckheere-Terpstra}: a trend test for ordinal variables. The
      asymptotic normal test (\code{jt}), the Monte Carlo permutation test
      (\code{mc-jt}) and the sequential Monte Carlo permutation test
      (\code{smc-jt}) are implemented.

  }

  Available conditional independence tests (and the respective labels) for
  \emph{Gaussian Bayesian networks} (normal variables) are:

  \itemize{

    \item \emph{linear correlation}: Pearson's linear correlation. The exact
      Student's t test (\code{cor}), the Monte Carlo permutation test
      (\code{mc-cor}) and the sequential Monte Carlo permutation test
      (\code{smc-cor}) are implemented.
    \item \emph{Fisher's Z}: a transformation of the linear correlation with
	  asymptotic normal distribution. The asymptotic normal test (\code{zf}),
	  the Monte Carlo permutation test (\code{mc-zf}) and the sequential Monte
	  Carlo permutation test (\code{smc-zf}) are implemented.
    \item \emph{mutual information}: an information-theoretic distance measure.
	  Again it is proportional to the log-likelihood ratio (they differ by a
	  \eqn{2n}{2n} factor). The asymptotic \eqn{\chi^2}{chi-square} test
	  (\code{mi-g}), the Monte Carlo permutation test (\code{mc-mi-g}) and the
	  sequential Monte Carlo permutation test (\code{smc-mi-g}) are
	  implemented.
    \item \emph{shrinkage estimator} for the \emph{mutual information}
      (\code{mi-g-sh}): an improved asymptotic \eqn{\chi^2}{chi-square} test
      based on the James-Stein estimator for the mutual information.

      Ledoit O, Wolf M (2003). "Improved Estimation of the Covariance Matrix
      of Stock Returns with an Application to Portfolio Selection".
      \emph{Journal of Empirical Finance}, \strong{10}:603--621.

  }

  Available conditional independence tests (and the respective labels) for
  \emph{hybrid Bayesian networks} (mixed discrete and normal variables) are:

  \itemize{

    \item \emph{mutual information}: an information-theoretic distance measure.
	  Again it is proportional to the log-likelihood ratio (they differ by a
	  \eqn{2n}{2n} factor). Only the asymptotic \eqn{\chi^2}{chi-square} test
	  (\code{mi-cg}) is implemented.

  }

}
\keyword{package}
