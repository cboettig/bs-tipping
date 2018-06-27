#' Get the Root of the Fish Model
#' 
#' Numerically finds the roots of the fish model given starting values (adults, [sun]fish, and juveniles) and parameter (harvest of adults, qE)
#' 
#' @param x a vector of length 3 with names A0 (adult bass), F0 (planktivorous fish), and J0 (juvenile bass)
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


#' Get the Root of the Fish Model in 2D
#' 
#' Numerically finds the roots of the fish model given starting values (juveniles and [sun]fish) and parameter (harvest of adults, qE)
#' 
#' @param x a vector of length 2 with names J0 (juveniles) and F0 (planktivorous fish)
#' @param pars vector of parameters to pass to model
#' @param maxiter maximum number of iterations to use when finding root
#' @param ... additional arguments to pass to \code{\link{fishStep2D}}
#' @return numeric vector with roots
#' @examples
#' getRoot2D(c(F0=1, J0=1), pars=c(qE=0.65))
#' getRoot2D(c(F0=1, J0=1000), pars=c(qE=0.65))
#' getRoot2D(c(F0=1, J0=1), pars=c(qE=0.5))
#' getRoot2D(c(F0=1, J0=1000), pars=c(qE=0.5))
#' @export
getRoot2D <- function(x, pars, maxiter=1E3, ...){
	# make function check x and pars for qE
	# can be convenient if you want to pass 1 object that has both state variables and parameter
	# such as when iterating throw rows of matrix with many combinations of states and parameters
	# only does this for qE, not other potential parameters
	if(missing(pars)){
		pars <- c(qE=unname(x["qE"]))
	}else{
		if(!"qE"%in%names(pars)){
			pars <- c(qE=unname(x["qE"]), pars)
		}else{
			pars <- c(pars)
		}
	}
	tryCatch({
		rootSolve::multiroot(f=fishStep2D, start=x[c("J0","F0")], parms=pars, maxiter=maxiter, ctol=1E-9, rtol=1E-9, atol=1E-9, ...)$root
	}, warning=function(w)c(J0=NA,F0=NA))
}



#' Get the Root of the Trophic Cascade Model from a Data Frame
#' 
#' Finds the roots of the trophic cascade model given starting values (juvenile bass, planktivorous fish) and parameter (harvest)
#' 
#' @param initialValues a \code{data.frame} with 3 columns named J0 (juvenile bass), F0 (planktivorous fish), and qE (harvest)
#' @param maxiter maximum number of iterations
#' @param ... additional arguments to pass to \code{\link{fishStep2D}}
#' @return a matrix with columns for initial values of J0 and F0, roots, and the parameter qE. Each row of output corresponds to a row of the input (though output reorded to ascending qE). Original columns for state variables (J0 and F0) will be renamed to "init.J0" and "init.F0". The equilibria that would be approached from these initial values, i.e. the 'roots', will take on the column names of "J0" and "F0". Thus, the columns J0 & F0 will not have the same values in the input as in the output (output are equilibrium values), unless the input states were already at equilibrium.
#' @examples
#' J0v <- seq(0, 1.5E3, length=2)
#' F0v <- 20
#' qEv <- seq(0.15, 1.5, length=10)
#' iV <- expand.grid(J0=J0v, F0=F0v, qE=qEv)
#' roots <- getRoot2D_df(iV)
#' @export
getRoot2D_df <- function(initialValues, maxiter=1E3, ...){
	rootGrid <- t(apply(initialValues, 1, getRoot2D, maxiter=maxiter, ...))
	rootGrid_qE <- data.matrix(data.frame(qE=initialValues[,"qE"], init=initialValues[,c("J0","F0")], rootGrid))
	rootGrid_qE <- rootGrid_qE[order(rootGrid_qE[,"qE"]),] # i should really remove this ...
	return(rootGrid_qE)
}
