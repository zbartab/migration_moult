\name{ini_gen}
\alias{ini.gen}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generating series of ini files}
\description{
  Faclitate to produce a series of ini files.
}
\usage{
ini.gen(a.nev = "baseline")
}
\arguments{
  \item{a.nev}{contains the name of the base ini (\file{.bini}) file}
}
\details{
  The function read a base ini (\file{.bini}) file and generates a series of
  ini files for the program \code{moult} based on this file. The
  \file{.bini} file is a text file containing a variable descriptor
  followed by a "\t" then a value or a series of values which are
  separeted by spaces. Then the function generates ini files for each
  possible combination for the values in the series. The first line of
  \file{.bini} file is special becaues it contains the basename of the
  generated files and whether the values for feather types should be
  equal ("TRUE") or not ("FALSE"). 
}
\author{Z. Barta}
\seealso{ \code{\link{rand.ini}}}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{moult}% __ONLY ONE__ keyword per line
