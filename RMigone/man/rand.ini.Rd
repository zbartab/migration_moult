\name{moult_ini}
\alias{rand.ini}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generates a random initialization file for moult}
\description{
  \code{rand.ini} generates a random initialization file for the moult
  program based on a sample initialization file.
}
\usage{
rand.ini(nev = "moult", n = 100, tart = 0.2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{nev}{the name of the sample ini file}
  \item{n}{number of random ini files}
  \item{tart}{\code{tart} gives the range (as proportion of the value in
    the sample file) within the random values vary in the random files.
}
\details{
  The sample file is a simple text file each line of which contains a
  variable descriptor and a value separated by a "\t". The generated
  values randomly vary uniformly between \code{x}-\code{x}*\code{tart}
  and \code{x}+\code{x}*\code{tart}.

  The name of the random file generated randomly.
}
\value{
  Used for its side effects
}
\author{Z. Barta}
\note{No care is taken to ensure unique random file names.}

 ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(nev="moult",n=100,tart=0.2){
  ini.start <- read.table(paste(nev,".ini",sep=""),col.names=c("var","value"));
  for(j in 1:n) {
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{moult}% __ONLY ONE__ keyword per line
