\name{whitelists-blacklists}
\alias{whitelists-blacklists}
\alias{whitelists and blacklists}
\title{Whitelists and blacklists in structure learning}
\description{

  How whitelists and blacklists are used in structure learning.

}
\section{Constraint-based Algorithms}{

  Constraint-based algorithms support arc whitelisting and blacklisting as
  follows:

  \itemize{

    \item blacklisted arcs are never present in the learned graph.
    \item arcs whitelisted in one direction only (i.e.
      \eqn{A \rightarrow B}{A -> B} is whitelisted but
      \eqn{B \rightarrow A}{B -> A} is not) have the respective reverse arcs
      blacklisted, and are always present in the learned graph.
    \item arcs whitelisted in both directions (i.e. both
      \eqn{A \rightarrow B}{A -> B} and \eqn{B \rightarrow A}{B -> A} are
	  whitelisted) are present in the learned graph, but their direction is set
	  by the learning algorithm.

  }

  Any arc whitelisted and blacklisted at the same time is assumed to be
  whitelisted, and is thus removed from the blacklist.

}
\section{Score-based Algorithms}{

  Score-based algorithms support arc whitelisting and blacklisting as follows:

  \itemize{

    \item blacklisted arcs are never present in the learned graph.
    \item arcs can only be whitelisted in a single direction, and are always
      present in the learned graph; it is not possible to whitelist arcs in both
      directions.

  }

}
\section{Hybrid Algorithms}{

  Hybrid algorithms use constraint-based (or pairwise mutual information)
  algorithms in the \emph{restrict phase} and score-based algorithms in the
  \emph{maximize phase}. Hence whitelists and blacklists are supported as
  follows:

  \itemize{

    \item whitelists and blacklists should be specified for the algorithm used
      in the restrict phase.
    \item if the whitelist contains any undirected arc, its consistent extension
      is used instead in the maximize phase.

  }

}
\section{Pairwise Mutual Information Algorithms}{

  In algorithms that learn undirected graphs, such as ARACNE and Chow-Liu, arcs
  are treated as being whitelisted or blacklisted in both directions even if
  only one direction is listed in the whitelist or blacklist. Again blacklisted
  arcs are never present in the learned graph and whitelisted arcs are
  guaranteed to be present in the learned graph.

}
\keyword{package}
