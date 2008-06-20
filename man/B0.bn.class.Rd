\name{bn class}
\alias{bn-class}
\title{

  The bn class structure

}
\description{

  The structure of an object of the \code{bn} S3 class.

}
\details{

  An object of class \code{bn} is a list containing at least the
  following components:

  \itemize{

    \item \code{learning}: a list containing some information about
       the results of the learning algorithm. It's never changed
       afterward.

    \itemize{

      \item \code{nodes}: a list. Each element is named after a node
        and contains the following elements:

      \itemize{

        \item \code{mb}: the Markov blanket of the node (a vector of
          character strings).

        \item \code{nbr}: the neighbourhood of the node (a vector of
          character strings).

      }

      \item \code{arcs}: the arcs of the Bayesian network (a two-column
        matrix, whose columns are labeled \code{from} and \code{to}).

      \item \code{whitelist}: a sanitized copy of the \code{whitelist}
        parameter (a two-column matrix, whose columns are labeled
        \code{from} and \code{to}).

      \item \code{blacklist}: a sanitized copy of the \code{blacklist}
        parameter (a two-column matrix, whose columns are labeled
        \code{from} and \code{to}).

      \item \code{test}: the label of the conditional independence test
        or used by the learning algorithm (a character string). The
        label of the network score is used for score-based algorithms.

      \item \code{ntests}: the number of conditional independence tests
        used in the learning (an integer value).

      \item \code{algo}: the label of the algorithm used in the learning
        process (a character string), or "random/generated" for randomly
        generated or empty networks.

      \item \code{args}: a list. The values of the parameters of
        either the conditional tests or the scores used in the learning
        process. Only the relevant ones are stored, so this may be
        an empty list.

        \itemize{

          \item \code{alpha}: the target nominal type I error rate (a
            numerical value) of the conditional independence tests.

          \item \code{iss}: a positive numerical value, the imaginary
            sample size used by the \code{bge} and \code{bde} scores.

          \item \code{phi}: a character string, either \code{heckerman}
            or \code{bottcher}; used by the \code{bge} score.

          \item \code{k}: a positive numerical value, the penalty per
            parameter used by the \code{aic} and \code{bic} scores.

        }

    }

    \item \code{nodes}: a list. Each element is named after a node
        and contains the following elements:

    \itemize{

      \item \code{mb}: the Markov blanket of the node (a vector of
        character strings).

      \item \code{nbr}: the neighbourhood of the node (a vector of
        character strings).

      \item \code{parents}: the parents of the node (a vector of
        character strings).

      \item \code{children}: the children of the node (a vector of
        character strings).

    }

    \item \code{arcs}: the arcs of the Bayesian network (a two-column
      matrix, whose columns are labeled \code{from} and \code{to}).

  }

}
\author{ Marco Scutari }
\keyword{ internal }
