\name{test counter}
\alias{test.counter}
\alias{increment.test.counter}
\alias{reset.test.counter}
\title{Manipulating the test counter}
\description{

  Check, increment or reset the test/score counter used in structure learning
  algorithms.

}
\usage{
test.counter()
increment.test.counter(i = 1)
reset.test.counter()
}
\arguments{
  \item{i}{a numeric value, which is added to the test counter.}
}
\value{

  A numeric value, the current value of the test counter.

}
\examples{
data(learning.test)
hc(learning.test)
test.counter()
reset.test.counter()
test.counter()
}
\author{ Marco Scutari }
\keyword{utilities}
