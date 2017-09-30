\name{plot_moult}
\alias{qp}
\alias{qp.sim}
\alias{qp.C}
\alias{qp.ser}
\alias{qpC.ser}
\alias{p.ser}
\title{Plotting moult policies}
\description{
  These functions are used to plotting the behaviour, policy and the
  reserves for moult runs.
}
\usage{
qp(nev = "moult", ee = 2, ll = 1)
qp.sim(y = 3, nev = "moult", ask = FALSE)
qp.C(y = 3, nev = "moult", ll = 1)
p.ser(y = 3, m.path = ".", patt = "\\.dat$",f.list=NULL)
qp.ser(m.path = ".", patt = "[^.][^C]\\.pop$",f.list=NULL)
qpC.ser(m.path = ".", patt = "\\.C\\.pop$",f.list=NULL)
}
\arguments{
  \item{nev}{The name of the file contains the data. Extension should
    not be given!
  \item{ee}{The number of "experience" class in the run}
  \item{ll}{The number of location in the run}
  \item{y}{The \code{y}th year should be shown}
  \item{ask}{Should the program wait after plotting?}
  \item{m.path}{The place where the data files are}
  \item{patt}{The matching pattern of the data file. The functions plot
    the data for each file matched against this pattern}
  \item{f.list}{contains the base names of the data files to plot. If it
    is NULL the \code{m.path} searched for files matching \code{patt}
    and those found plotted}
}
\details{
  \code{qp} show the behaviour based on the \code{OutPut} function of
  \code{moult}. \code{qp.sim} plots the behaviour and the reserves based
  on the simulation results. \code{qp.C} gives the behaviour for a
  cohort. \code{p.ser}, \code{qp.ser} and \code{qpC.ser} are wrappers
  for \code{qp.sim}, \code{qp} and \code{qp.C}, respectively, allowing the
  plotting for a series of files. 
}
\author{Z. Barta}
\keyword{utilities}% at least one, from doc/KEYWORDS
\keyword{moult}% __ONLY ONE__ keyword per line
\keyword{plotting}% __ONLY ONE__ keyword per line
