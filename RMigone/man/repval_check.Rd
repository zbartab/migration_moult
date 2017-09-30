\name{repval_check}
\alias{irott}
\alias{olvasott}
\title{Checking the correctness of repval array}
\description{
  These functions check whether the reproductive values in the
  repval array of moult program is written/read before they are read/written
}
\usage{
irott(i, x, data)
olvasott(i, x, data)
}
\arguments{
  \item{i}{\code{i} is ID of the repval to check here}
  \item{x}{\code{x} is the data.frame containing all the repval to check}
  \item{data}{\code{data} is the data.frame containing the values to
    which compare \code{x}}
}
\details{
  Both function search the \code{data} for records which match with
  \code{x[i,]}. If at least a match found then prints \code{i}, the ID
  of \code{x}, the ID of matched record in \code{data} and how many
  match occured.  
}
\value{
  Returns the number of matched records.
}
\author{Z. Barta}
\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }
\keyword{utilities}
\keyword{moult}
