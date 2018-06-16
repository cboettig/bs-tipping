#' Get the Root of the Fish Model
#' Numerically finds the roots of the fish model given starting values (adults, [sun]fish, and juveniles) and parameter (harvest of adults, qE)
#' @param x a vector of length 3 with names A0 (water P), F0 (mud P), and J0 (P loading)
#' @param pars vector of parameters to pass to model
#' @param maxiter maximum number of iterations to use when finding root
#' @param ... additional arguments to pass to \code{\link{fishStep}}
#' @return numeric vector with roots
#' @examples
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.65))
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.65))
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.5))
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.5))
#' @export
getRoot <- function(x, pars, maxiter=1E3, ...){
	# make function check x and pars for qE
	# can be convenient if you want to pass 1 object that has both state variables and parameter
	# such as when iterating throw rows of matrix with many combinations of states and parameters
	# only does this for qE, not other potential parameters
	if(missing(pars) | is.null(pars)){
		pars <- c(qE=unname(x["qE"]))
	}else{
		if(!"qE"%in%names(pars)){
			pars <- c(qE=unname(x["qE"]), pars)
		}else{
			pars <- c(pars)
		}
	}
	tryCatch({
		rootSolve::multiroot(f=fishStep, start=x[c("A0","F0","J0")], parms=pars, maxiter=maxiter, ctol=1E-9, rtol=1E-9, atol=1E-9, ...)$root
	}, warning=function(w)c(A0=NA,F0=NA,J0=NA))
}