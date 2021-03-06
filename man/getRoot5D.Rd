% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRoot.R
\name{getRoot5D}
\alias{getRoot5D}
\title{Get the Root of the Full Ecosystem Model in 5D}
\usage{
getRoot5D(x, pars, maxiter = 1000, ...)
}
\arguments{
\item{x}{a vector of length 5 with names A0 (adult bass), F0 (planktivorous fish),  J0 (juvenile bass), H0 (herbivorous zooplankton), and P0 (phytoplankton)}

\item{pars}{vector of parameters to pass to model}

\item{maxiter}{maximum number of iterations to use when finding root}

\item{...}{additional arguments to pass to \code{\link{ecoStep}}}
}
\value{
numeric vector with roots
}
\description{
Numerically finds the roots of the ecosystem model given starting values (adults, [sun]fish, juveniles, herbivores, and phytoplankton) and parameter (harvest of adults, qE)
}
\examples{
getRoot5D(c(A0=15, F0=50, J0=10, H0=100, P0=100), pars=c(qE=-0.02)) # stocking, high starting A
getRoot5D(c(A0=2, F0=50, J0=5, H0=100, P0=100), pars=c(qE=-0.02)) # stocking, but start low A
getRoot5D(c(A0=15, F0=1, J0=15, H0=100, P0=100), pars=c(qE=0.02)) # harvesting, start high A
}

