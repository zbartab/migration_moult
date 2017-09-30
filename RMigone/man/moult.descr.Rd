\name{moult.descr}
\alias{moult.descr}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generates some descriptive statistics}
\description{
  Take a data frame containing the results for \emph{one} year of a
  \code{moult} run and calculate some statistics from it.
}
\usage{
moult.descr(ds)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ds}{the data frame to contains \code{moult} simulation results
    of \emph{one} year.} 
}
\details{
  The data frame is created as \code{ds <-
    read.table("moult.dat",header=TRUE)}. 
}
\value{
  Returns the calculated statistics as a vector:
  \item{1}{mean number of broods}
  \item{2}{mean number of moulting feather 1}
  \item{3}{mean number of moulting feather 2}
  \item{4}{mean length (in weeks) of moulting feather 1}
  \item{5}{mean length (in weeks) of moulting feather 2}
  \item{6}{mean length of moult}
  \item{7}{mean length of breeding}
  \item{8}{mean timing of the first brood}
}
\author{Z. Barta}
\note{Under development}
\seealso{ \code{\link{qp}}, \code{\link{qp.sim}}, \code{\link{qp.C}} }

\keyword{statistics}% at least one, from doc/KEYWORDS
\keyword{moult}% __ONLY ONE__ keyword per line
