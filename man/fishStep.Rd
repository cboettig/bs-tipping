% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fishStep.R
\name{fishStep}
\alias{fishStep}
\title{Fish Step}
\usage{
fishStep(X, pars = c(qE = 1))
}
\arguments{
\item{X}{vector of length 3 with names A0, F0, and J0.  Indicates starting value of state variables}

\item{pars}{a named vector with at least qE, indicating harvesting rate of adults; can include other relevant parameters, see \code{\link{ecoStep}}}
}
\value{
named vector of length 3, with names corresponding to juvenile bass (J0), adult bass (A0), and sunfish (F0) abundances.
}
\description{
Compute the rate of change of the three fish state variables given initial conditions and parameters
}
\examples{
fishStep(X=c(A0=1.5,F0=25,J0=2), pars=c(qE=0.1))
}

