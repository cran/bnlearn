\name{network-scores}
\alias{network-scores}
\alias{network scores}
\title{Network scores}
\description{

  Overview of the network scores implemented in \pkg{bnlearn}, with the
  respective reference publications.

}
\details{

  Note that all instances of AIC, BIC and eBIC are rescaled by \code{-2} to make
  them comparable with the corresponding marginal log-likelihoods. The same is
  true for the penalised node average log-likelihoods.

  Available scores (and the respective labels) for \emph{discrete Bayesian
  networks} (categorical variables) are:

  \itemize{

    \item The multinomial \emph{log-likelihood} (\code{"loglik"}) score, which is
      equivalent to the \emph{entropy measure} used in Weka.
    \item The \emph{Akaike Information Criterion} (AIC) score (\code{"aic"}).
    \item The \emph{Bayesian Information Criterion} (BIC) score (\code{"bic"}),
      which is equivalent to the \emph{Minimum Description Length} (MDL) and is
      also known as \emph{Schwarz Information Criterion}.

      Chickering DM (1995). "A Transformational Characterization of Equivalent
        Bayesian Network Structures." \emph{Proceedings of the Eleventh Annual
        Conference on Uncertainty in Artificial Intelligence}, 87--98.

    \item The \emph{extended Bayesian Information Criterion} (\code{"ebic"}),
      which adds a second penalty to BIC to penalise dense networks.

      Foygel R, Drton M (2010). "Extended Bayesian Information Criteria for
        Gaussian Graphical Models." NIPS 23, 604--612.

    \item The \emph{predictive log-likelihood} (\code{"pred-loglik"}) computed on
      a separate test set.

      Chickering DM, Heckerman D (2000). "A Comparison of Scientific and
        Engineering Criteria for Bayesian Model Selection." \emph{Statistics
        and Computing}, \strong{10}:55--62.

      Scutari M, Vitolo C, Tucker A (2019). "Learning Bayesian Networks from
        Big Data with Greedy Search: Computational Complexity and Efficient
        Implementation." \emph{Statistics and Computing},
        \strong{25}(9):1095--1108.

    \item The logarithm of the \emph{Bayesian Dirichlet equivalent (uniform)}
      score (\code{"bde"}) (also denoted BDeu), a score-equivalent Dirichlet
      posterior density.

      Heckerman D, Geiger D, Chickering DM (1995). "Learning Bayesian Networks:
        The Combination of Knowledge and Statistical Data."  \emph{Machine
        Learning}, \strong{20}(3):197--243.

      Castelo R, Siebes A (2000). "Priors on Network Structures. Biasing the
        Search for Bayesian Networks." \emph{International Journal of
        Approximate	Reasoning}, \strong{24}(1):39--57.

    \item The logarithm of the \emph{Bayesian Dirichlet sparse} score
	  (\code{"bds"}) (BDs), a sparsity-inducing Dirichlet posterior density (not
	  score equivalent).

      Scutari M (2016). "An Empirical-Bayes Score for Discrete Bayesian
        Networks." \emph{Journal of Machine Learning Research},
        \strong{52}:438--448.

    \item The logarithm of the \emph{Bayesian Dirichlet} score with
      \emph{Jeffrey's prior} (not score equivalent).

      Suzuki J (2016). "A Theoretical Analysis of the BDeu Scores in Bayesian
        Network Structure Learning." \emph{Behaviormetrika},
        \strong{44}(1):97--116.

    \item The logarithm of the modified \emph{Bayesian Dirichlet equivalent}
	  score (\code{"mbde"}) for mixtures of experimental and observational data
	  (not score equivalent).

      Cooper GF, Yoo C (1999). "Causal Discovery from a Mixture of Experimental
        and Observational Data." \emph{Proceedings of the Fifteenth Annual
        Conference on Uncertainty in Artificial Intelligence}, 116--125.

    \item The logarithm of the \emph{K2} score (\code{"k2"}), a Dirichlet
      posterior density (not score equivalent).

      Korb K, Nicholson AE (2010). \emph{Bayesian Artificial Intelligence}.
        Chapman & Hall/CRC, 2nd edition.

    \item The logarithm of the \emph{factorized normalized maximum likelihood}
      score (\code{"fnml"}, not score equivalent).

      Silander T, Roos T, Kontkanen P, Myllymaki P (2008). "Factorized
        Normalized Maximum Likelihood Criterion for Learning Bayesian Network
        Structures." \emph{Proceedings of the 4th European Workshop on
        Probabilistic Graphical Models}, 257--272.

    \item The logarithm of the \emph{quotient normalized maximum likelihood}
      (\code{"qnml"}).

      Silander T, Leppa-Abo J, Jaasaari, Roos T (2018). "Quotient Normalized
        Maximum Likelihood Criterion for Learning Bayesian Network Structures."
        \emph{Proceedings of Machine Learning Research}, \strong{84}:948--957.

    \item The node-average (log-)likelihood (\code{"nal"}) and the penalised
      node-average (log-)likelihood (\code{"pnal"}).

      Bodewes T, Scutari M (2021). "Learning Bayesian Networks from Incomplete
        Data with the Node-Averaged Likelihood." \emph{International Journal of
        Approximate Reasoning}, \strong{138}:145--160.

  }

  Available scores (and the respective labels) for \emph{Gaussian Bayesian
  networks} (normal variables) are:

  \itemize{

    \item The multivariate Gaussian \emph{log-likelihood} (\code{"loglik-g"})
      score.
    \item The corresponding \emph{Akaike Information Criterion} (AIC) score
      (\code{"aic-g"}).
    \item The corresponding \emph{Bayesian Information Criterion} (BIC) score
      (\code{"bic-g"}).

      Geiger D, Heckerman D (1994). "Learning Gaussian Networks."
        \emph{Proceedings of the Tenth Annual Conference on Uncertainty in
        Artificial Intelligence}, 235--243.

    \item The \emph{extended Bayesian Information Criterion} (\code{"ebic-g"}),
      which adds a second penalty to BIC to penalise dense networks. Like BIC,
      eBIC is rescaled by \code{-2}.

      Foygel R, Drton M (2010). "Extended Bayesian Information Criteria for
        Gaussian Graphical Models." NIPS 23, 604--612.

    \item The \emph{predictive log-likelihood} (\code{"pred-loglik-g"}) computed
      on a separate test set. The reference paper is the same as that for
      \code{"pred-loglik"}. It is currently implemented to be score-equivalent
      like \code{"pred-loglik"}, but this may change.

    \item A score-equivalent \emph{Gaussian posterior density} (\code{"bge"}).

      Kuipers J, Moffa G, Heckerman D (2014). "Addendum on the Scoring of
        Gaussian Directed Acyclic Graphical Models." \emph{The Annals of
        Statistics}, \strong{42}(4):1689--1691.

    \item The node-average (log-)likelihood (\code{"nal-g"}) and the penalised
      node-average (log-)likelihood (\code{"pnal-g"}).

      Bodewes T, Scutari M (2021). "Learning Bayesian Networks from Incomplete
        Data with the Node-Averaged Likelihood." \emph{International Journal of
        Approximate Reasoning}, \strong{138}:145--160.

  }

  Available scores (and the respective labels) for \emph{conditional Gaussian
  Bayesian networks} (mixed categorical and normal variables) are:

  \itemize{

    \item The conditional linear Gaussian \emph{log-likelihood}
      (\code{"loglik-cg"}) score.
    \item The corresponding \emph{Akaike Information Criterion} (AIC) score
      (\code{"aic-cg"}).
    \item The corresponding \emph{Bayesian Information Criterion} (BIC) score
      (\code{"bic-cg"}).

    \item The \emph{extended Bayesian Information Criterion} (\code{"ebic-cg"}),
      which adds a second penalty to BIC to penalise dense networks.

      Foygel R, Drton M (2010). "Extended Bayesian Information Criteria for
        Gaussian Graphical Models." NIPS 23, 604--612.

    \item The \emph{predictive log-likelihood} (\code{"pred-loglik-cg"}) computed
      on a separate test set. The reference paper is the same as that for
      \code{"pred-loglik"}.

    \item The node-average (log-)likelihood (\code{"nal-cg"}) and the penalised
      node-average (log-)likelihood (\code{"pnal-cg"}).

      Bodewes T, Scutari M (2021). "Learning Bayesian Networks from Incomplete
        Data with the Node-Averaged Likelihood." \emph{International Journal of
        Approximate Reasoning}, \strong{138}:145--160.

  }

  Available scores (and the respective labels) for \emph{zero-inflated Bayesian
  networks} (non-integer count variables with zero inflation) are:

  \itemize{

    \item The zero-inflated \emph{hyper-Poisson log-likelihood}
      (\code{"loglik-zihp"}) and \emph{negative binomial log-likelihood}
      (\code{"loglik-zinb"}) scores.
    \item The corresponding \emph{Akaike Information Criterion} (AIC) scores
      (\code{"aic-zihp"} and \code{"aic-zinb"}).
    \item The corresponding \emph{Bayesian Information Criterion} (BIC) scores
      (\code{"bic-zihp"} and \code{"bic-zinb"}).

      Choi J, Ni Y (2023). "Model-Based Causal Discovery for Zero-Inflated Count
        Data." \emph{Journal of Machine Learning Research},
        \strong{24}(200):1--32.

    \item The corresponding node-average (log-)likelihood (\code{"nal-zihp"} and
      \code{"nal-zinb"}) and the penalised node-average (log-)likelihood
      (\code{"pnal-zihp"} and \code{"pnal-zinb"}) scores.

      Bodewes T, Scutari M (2021). "Learning Bayesian Networks from Incomplete
        Data with the Node-Averaged Likelihood." \emph{International Journal of
        Approximate Reasoning}, \strong{138}:145--160.

    The local distributions underlying these scores are estimated by
    expectation-maximization with the generalised linear model components fitted
    by iteratively reweighted least squares (as in \code{\link{bn.fit}}); the
    estimator is controlled by the \code{em.max.iter}, \code{em.tol} and
    \code{m.step} arguments documented (with their defaults) in
    \code{\link{score}}.

  }

  Other scores (and the respective labels):

  \itemize{

    \item A \emph{custom decomposable} (\code{custom}) score interface that
      takes an R function as an argument. It can be used to trial experimental
      score functions without having to code them in C and hook them up to the
      internals of \pkg{bnlearn}.

  }

}
\keyword{package}
