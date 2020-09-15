\name{network-scores}
\alias{network-scores}
\alias{network scores}
\title{Network scores}
\description{

  Overview of the network scores implemented in \pkg{bnlearn}, with the
  respective reference publications.

}
\details{

  Available scores (and the respective labels) for \emph{discrete Bayesian
  networks} (categorical variables) are:

  \itemize{

    \item the multinomial \emph{log-likelihood} (\code{loglik}) score, which is
      equivalent to the \emph{entropy measure} used in Weka.
    \item the \emph{Akaike Information Criterion} score (\code{aic}).
    \item the \emph{Bayesian Information Criterion} score (\code{bic}), which is
      equivalent to the \emph{Minimum Description Length} (MDL) and is also
      known as \emph{Schwarz Information Criterion}.

      Chickering DM (1995). "A Transformational Characterization of Equivalent
		Bayesian Network Structures". \emph{Proceedings of the Eleventh Annual
		Conference on Uncertainty in Artificial Intelligence}, 87--98.

    \item the \emph{predictive log-likelihood} (\code{pred-loglik}) computed on
      a separate test set.

      Chickering DM, Heckerman D (2000). "A Comparison of Scientific and
		Engineering Criteria for Bayesian Model Selection". \emph{Statistics
		and Computing}, \strong{10}:55--62.

      Scutari M, Vitolo C, Tucker A (2019). "Learning Bayesian Networks from
		Big Data with Greedy Search: Computational Complexity and Efficient
		Implementation". \emph{Statistics and Computing}, online first.

    \item the logarithm of the \emph{Bayesian Dirichlet equivalent (uniform)}
      score (\code{bde}) (also denoted BDeu), a score equivalent Dirichlet
      posterior density.

      Heckerman D, Geiger D, Chickering DM (1995). "Learning Bayesian Networks:
		The Combination of Knowledge and Statistical Data".  \emph{Machine
		Learning}, \strong{20}(3):197--243.

      Castelo R, Siebes A (2000). "Priors on Network Structures. Biasing the
		Search for Bayesian Networks". \emph{International Journal of
		Approximate	Reasoning}, \strong{24}(1):39--57.

    \item the logarithm of the \emph{Bayesian Dirichlet sparse} score
	  (\code{bds}) (BDs), a sparsity-inducing Dirichlet posterior density (not
	  score equivalent).

      Scutari M (2016). "An Empirical-Bayes Score for Discrete Bayesian
		Networks". \emph{Journal of Machine Learning Research},
		\strong{52}:438--448.

    \item the logarithm of the \emph{Bayesian Dirichlet} score with
      \emph{Jeffrey's prior} (not score equivalent).

      Suzuki J (2016). "A Theoretical Analysis of the BDeu Scores in Bayesian
		Network Structure Learning". \emph{Behaviormetrika},
		\strong{44}(1):97--116.

    \item the logarithm of the modified \emph{Bayesian Dirichlet equivalent}
	  score (\code{mbde}) for mixtures of experimental and observational data
	  (not score equivalent).

      Cooper GF, Yoo C (1999). "Causal Discovery from a Mixture of Experimental
		and Observational Data". \emph{Proceedings of the Fifteenth Annual
		Conference on Uncertainty in Artificial Intelligence}, 116--125.

    \item the logarithm of the \emph{locally averaged Bayesian Dirichlet} score
      (\code{bdla}, not score equivalent).

      Cano A, Gomez-Olmedo M, Masegosa AR, Moral S (2013). "Locally Averaged
		Bayesian Dirichlet Metrics for Learning the Structure and the
		Parameters of Bayesian Networks". \emph{International Journal of
		Approximate Reasoning}, \strong{54}:526--540.

    \item the logarithm of the \emph{K2} score (\code{k2}), a Dirichlet
      posterior density (not score equivalent).

      Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
        Chapman & Hall/CRC, 2nd edition.

  }

  Available scores (and the respective labels) for \emph{Gaussian Bayesian
  networks} (normal variables) are:

  \itemize{

    \item the multivariate Gaussian \emph{log-likelihood} (\code{loglik-g})
      score.
    \item the corresponding \emph{Akaike Information Criterion} score
      (\code{aic-g}).
    \item the corresponding \emph{Bayesian Information Criterion} score
      (\code{bic-g}).

      Geiger D, Heckerman D (1994). "Learning Gaussian Networks".
        \emph{Proceedings of the Tenth Annual Conference on Uncertainty in
        Artificial Intelligence}, 235--243.

    \item the \emph{predictive log-likelihood} (\code{pred-loglik-g}) computed
	  on a separate test set. The reference paper is the same as that for
	  \code{pred-loglik}. It is currently implemented to be score-equivalent
      like \code{pred-loglik}, but that may be subject to change.

    \item a score equivalent \emph{Gaussian posterior density} (\code{bge}).

      Kuipers J, Moffa G, Heckerman D (2014). "Addendum on the Scoring of
        Gaussian Directed Acyclic Graphical Models". \emph{The Annals of
        Statistics}, \strong{42}(4):1689--1691.

  }

  Available scores (and the respective labels) for \emph{hybrid Bayesian
  networks} (mixed categorical and normal variables) are:

  \itemize{

    \item the conditional linear Gaussian \emph{log-likelihood}
      (\code{loglik-cg}) score.
    \item the corresponding \emph{Akaike Information Criterion} score
    \code{aic-cg}).
    \item the corresponding \emph{Bayesian Information Criterion} score
      (\code{bic-cg}).
    \item the \emph{predictive log-likelihood} (\code{pred-loglik-cg}) computed
      on a separate test set. The reference paper is the same as that for
      \code{pred-loglik}.

  }

  Other scores (and the respective labels):

  \itemize{

    \item a \emph{custom decomposable} (\code{custom}) score interface that
	  takes an R function as an argument. It can be used to trial experimental
	  score functions without having to code them in C and hook them up to the
	  internals of \pkg{bnlearn}.

  }

}
\keyword{package}
