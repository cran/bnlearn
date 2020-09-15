\name{structure-learning}
\alias{structure-learning}
\alias{structure learning}
\title{Structure learning algorithms}
\description{

  Overview of the structure learning algorithms implemented in \pkg{bnlearn},
  with the respective reference publications.

}
\section{Available Constraint-Based Learning Algorithms}{

  \itemize{

    \item \emph{PC} (\code{\link{pc.stable}}), a modern implementation of the
      first practical constraint-based structure learning algorithm.

      Colombo D, Maathuis MH (2014). "Order-Independent Constraint-Based Causal
        Structure Learning". \emph{Journal of Machine Learning Research},
        \strong{15}:3921--3962.

    \item \emph{Grow-Shrink} (\code{\link{gs}}): based on the \emph{Grow-Shrink
      Markov Blanket}, the first (and simplest) Markov blanket detection
      algorithm used in a structure learning algorithm.

      Margaritis D (2003). \emph{Learning Bayesian Network Model Structure from
        Data}. Ph.D. thesis, School of Computer Science, Carnegie-Mellon
        University, Pittsburgh, PA.

    \item \emph{Incremental Association} (\code{\link{iamb}}): based on the
      Markov blanket detection algorithm of the same name, which is based on
      a two-phase selection scheme (a forward selection followed by an attempt
      to remove false positives).

      Tsamardinos I, Aliferis CF, Statnikov A (2003). "Algorithms for Large
		Scale Markov Blanket Discovery". \emph{Proceedings of the Sixteenth
		International Florida Artificial Intelligence Research Society
		Conference}, 376--381.

    \item \emph{Fast Incremental Association} (\code{\link{fast.iamb}}): a
      variant of IAMB which uses speculative stepwise forward selection to
      reduce the number of conditional independence tests.
    \item \emph{Interleaved Incremental Association} (\code{\link{inter.iamb}}):
      another variant of IAMB which uses forward stepwise selection to avoid
      false positives in the Markov blanket detection phase.

      Yaramakala S, Margaritis D (2005). "Speculative Markov Blanket Discovery
		for Optimal Feature Selection". \emph{Proceedings of the Fifth IEEE
		International Conference on Data Mining}, 809--812.

    \item \emph{Incremental Association with FDR} (\code{\link{iamb.fdr}}): a
      variant of IAMB which adjusts the tests significance threshold with FDR.

      Pena JM (2008). "Learning Gaussian Graphical Models of Gene Networks with
        False Discovery Rate Control". \emph{Proceedings of the Sixth European
        Conference on Evolutionary Computation, Machine Learning and Data Mining
        in Bioinformatics}, 165--176.

      Gasse M, Aussem A, Elghazel H (2014). "A Hybrid Algorithm for Bayesian
        Network Structure Learning with Application to Multi-Label Learning".
        \emph{Expert Systems with Applications}, \strong{41}(15):6755--6772.

  }

  \pkg{bnlearn} includes two implementations of each algorithm: a vanilla
  implementation, and a parallel one that requires a running cluster set
  up with the \code{makeCluster} function from the \pkg{parallel} package.

}
\section{Available Score-based Learning Algorithms}{

  \itemize{

    \item \emph{Hill-Climbing} (\code{\link{hc}}): a \emph{hill climbing}
	  greedy search that explores the space of the directed acyclic graphs by
	  single-arc addition, removal and reversals; with random restarts to avoid
	  local optima. The optimized implementation uses score caching, score
	  decomposability and score equivalence to reduce the number of duplicated
	  tests.
    \item \emph{Tabu Search} (\code{\link{tabu}}): a modified hill-climbing
      able to escape local optima by selecting a network that minimally
      decreases the score function.

      Russell SJ, Norvig P (2009). \emph{Artificial Intelligence: A Modern
        Approach}. Prentice Hall, 3rd edition.

  }

}
\section{Available Hybrid Learning Algorithms}{

  \itemize{

    \item \emph{Max-Min Hill-Climbing} (\code{\link{mmhc}}): a hybrid algorithm
      which combines the Max-Min Parents and Children algorithm (to restrict the
      search space) and the Hill-Climbing algorithm (to find the optimal network
      structure in the restricted space).

      Tsamardinos I, Brown LE, Aliferis CF (2006). "The Max-Min Hill-Climbing
        Bayesian Network Structure Learning Algorithm". \emph{Machine Learning},
        \strong{65}(1):31--78.

    \item \emph{Restricted Maximization} (\code{\link{rsmax2}}): a general
      implementation of the Sparse Candidate algorithms, which can use any
      combination of constraint-based and score-based algorithms.

      Friedman N, Nachman I, Pe'er D (1999). "Learning Bayesian Network
        Structure from Massive Datasets: the Sparse Candidate Algorithm."
        \emph{Proceedings of the Fifteenth Conference on Uncertainty in
        Artificial Intelligence (UAI)}, 206--215.

    \item \emph{Hybrid HPC} (\code{\link{h2pc}}): a hybrid algorithm combining
      HPC and hill-climbing.

      Gasse M, Aussem A, Elghazel H (2014). "A Hybrid Algorithm for Bayesian
		Network Structure Learning with Application to Multi-Label Learning".
		\emph{Expert Systems with Applications}, \strong{41}(15):6755--6772.

  }

}
\section{Other (Constraint-Based) Local Discovery Algorithms}{

  These algorithms learn the structure of the undirected graph underlying the
  Bayesian network, which is known as the \emph{skeleton} of the network.
  Therefore by default all arcs are undirected, and no attempt is made to
  detect their orientation. They are often used in hybrid learning algorithms.

  \itemize{

    \item \emph{Max-Min Parents and Children} (\code{\link{mmpc}}): a forward
      selection technique for neighbourhood detection based on the maximization
      of the minimum association measure observed with any subset of the nodes
      selected in the previous iterations.

      Tsamardinos I, Aliferis CF, Statnikov A (2003). "Time and Sample Efficient
		Discovery of Markov Blankets and Direct Causal Relations".
		\emph{Proceedings of the Ninth ACM SIGKDD International Conference on
		Knowledge Discovery and Data Mining}, 673--678.

    \item \emph{Hiton Parents and Children} (\code{\link{si.hiton.pc}}): a fast
      forward selection technique for neighbourhood detection designed to
      exclude nodes early based on the marginal association. The implementation
      follows the Semi-Interleaved variant of the algorithm.

      Aliferis FC, Statnikov A, Tsamardinos I, Subramani M, Koutsoukos XD (2010).
		"Local Causal and Markov Blanket Induction for Causal Discovery and
		Feature Selection for Classification Part I: Algorithms and Empirical
		Evaluation".  \emph{Journal of Machine Learning Research},
		\strong{11}:171--234.

    \item \emph{Hybrid Parents and Children} (\code{\link{hpc}}): an algorithm
      building on \code{iamb.fdr} to learn the parents and children of each node
      like \code{mmpc} and \code{si.hiton.pc}. The reference publication is the
      same as that for Hybrid HPC.

  }

}
\section{Pairwise Mutual Information Algorithms}{

  These algorithms learn approximate network structures using only pairwise
  mutual information.

  \itemize{

    \item \emph{Chow-Liu} (\code{\link{chow.liu}}): an application of the
      minimum-weight spanning tree and the information inequality. It learns
      the tree structure closest to the true one in the probability space.

      Chow CK, Liu CN (1968). "Approximating Discrete Probability Distributions
        with Dependence Trees", IEEE Transactions on Information Theory, IT-14
        \strong{3}:462--467.

    \item \emph{ARACNE} (\code{\link{aracne}}): an improved version of the
      Chow-Liu algorithm that is able to learn polytrees.

      Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R,
		Califano A (2006). "ARACNE: An Algorithm for the Reconstruction of Gene
		Regulatory Networks in a Mammalian Cellular Context". \emph{BMC
		Bioinformatics}, \strong{7}(Suppl 1):S7.

  }

  All these algorithms have two implementations (vanilla and parallel) like
  other constraint-based algorithms.

}
\keyword{package}
