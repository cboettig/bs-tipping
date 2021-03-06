% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability.R
\name{stateMotion}
\alias{stateMotion}
\title{Equations of Motion}
\usage{
stateMotion(state0, stateVar = c("A0", "F0", "J0"), pars)
}
\arguments{
\item{state0}{numeric, value of a state variable (can be scalar or vector, not named)}

\item{stateVar}{character, name of state variable whose motion you want}

\item{pars}{named vector of parameters, needs to at least include qE}
}
\value{
numeric scalar or vector indicating the rate of change of the state variable
}
\description{
Provides the motion of each state variable (change in that variable as a function of itself and parameters)
}
\note{
These functions should essentially replace \code{\link{dFJ_dt_1state}}. The A motion function was not an option in \code{\link{dFJ_dt_1state}}, though. Both A motion and J motion are done without much algebra or substitution; however, getting the equation/ function for F motion requires doing some alebraic substitution to get the A state variable to cancel out. Dividing by A and letting it cancel probably involves some extra assumption that might make the F motion equation/ function invalid in certain areas of state-space.
}

