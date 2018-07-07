
#' Calculate the Roots of F
#' 
#' Calculate the roots of F (planktivorous fish) by solving a polynomial
#' 
#' @param pars named vector of parameters found in \code{\link{ecoStep}}. At a minimum, qE needs to be supplied. Optionally, other parameter arguments can be supplied here, but if they aren't, they are taken from the formals (defaults) of \code{\link{ecoStep}}.
#' 
#' @return numeric value(s) corresponding to roots of F
#' @examples
#' calcRoot_F(pars=c(qE=0.5))
calcRoot_F <- function(pars){
	parsF <- unlist(formals(ecoStep)[c("fA", "cJA", "cJF", "cFA", "vuln", "hide", "surv", "Fo", "DF")])
	
	if(missing(pars)){ # if a function requires qE, pars needs to be supplied (e.g., pars=c(qE=0.1))
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	
	with(as.list(pars), {
		qE <- pars['qE']
		a <- cFA
		x <- cJA
		z <- cJF
		f <- fA
		v <- vuln
		h <- hide
		s <- surv
		d <- DF
	
		# polyCoef0 <- -d*Fo
		# polyCoef1 <- (v+h)*(d*x*(qE+1-s)+a*s*f)
		# polyCoef2 <- z*((d*x-a*v)*(qE+1-s)+a*s*f)
		
		polyCoef0 <- 1
		polyCoef1 <- 1
		polyCoef2 <- 1
	
		polyroot(c(polyCoef0, polyCoef1, polyCoef2))
	})
	
	
}



#' Temporal dynamics of F or J
#' 
#' Temporal dynamics of F or J, but as a function of F or J and parameters (not other state variables)
#' 
#' @param State0 abundance of planktivorous fish (F0) or abundance of juvenile bass (J0)
#' @param pars named vector of parameters found in \code{\link{ecoStep}}. At a minimum, qE needs to be supplied. Optionally, other parameter arguments can be supplied here, but if they aren't, they are taken from the formals (defaults) of \code{\link{ecoStep}}.
#' @param stateName a length-1 character indicating whether the value supplied (and value returned) pertains to "F0" (planktivores) or "J0" (juvenile bass). Juveniles ("J0") is chosen by default.
#' @param parts Logical, if FALSE (default), function returns the rate of change in the indicated state variable. If TRUE, function returns the rate of 'growth' (>=0) and the rate of 'consumption' (<=0) for the indicated state variable.
#' 
#' @return numeric value indicating temporal rate of change of F or J; or, if parts==TRUE, a named vector of length 2 indicate the growth and consumption components of the temporal rate of change of F or J.
#' @examples
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.5)) # find equilibria when bass are rare
#' dFJ_dt_1state(State0=17.890184, pars=c(qE=0.5), stateName="F0") # F0 set to equilibrium when bass start as rare
#' dFJ_dt_1state(State0=3.06, pars=c(qE=0.5), stateName="J0") # J0 set to equilibrium when bass start as rare
#' 
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.5)) # find equilibria when bass are abundant
#' dFJ_dt_1state(State0=0.06712064, pars=c(qE=0.5), stateName="F0") # F0 set to equilibrium when bass start as abundant
#' dFJ_dt_1state(State0=992.5699, pars=c(qE=0.5), stateName="J0") # J0 set to equilibrium when bass start as abundant
#' 
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.65)) # find equilibria when bass are rare
#' dFJ_dt_1state(State0=100, pars=c(qE=0.65), stateName="F0") # F0 set to equilibrium when bass start as rare # WTF?!
#' fishStep(X=c(A0=6.275129e-18, F0=100, J0=1.443280e-17))
#' dFJ_dt_1state(State0=1.44e-17, pars=c(qE=0.65), stateName="J0") # J0 set to equilibrium when bass start as rare # Makes sense
#' 
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.65)) # find equilibria when bass are abundant
#' dFJ_dt_1state(State0=0.09136212, pars=c(qE=0.65), stateName="F0") # check planktivore equation for rate of change
#' dFJ_dt_1state(State0=838.3849, pars=c(qE=0.65), stateName="J0") # check juvenile equation for rate of change
#' fishStep(X=c(A0=364.51517642, F0=0.09136212, J0=838.38490576)) # re-examine rates of change of fish near equilibrium
#' @export
dFJ_dt_1state <- function(State0, pars, stateName=c("J0","F0"), parts=FALSE){
	stateName <- match.arg(stateName)
	if(stateName=="J0"){
		J0 <- State0
	}
	if(stateName=="F0"){
		F0 <- State0
	}
	parsF <- unlist(formals(ecoStep)[c("fA", "cJA", "cJF", "cFA", "vuln", "hide", "surv", "Fo", "DF")])
	
	if(missing(pars)){ # if a function requires qE, pars needs to be supplied (e.g., pars=c(qE=0.1))
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	
	out <- with(as.list(pars), {
		qE <- pars['qE']
		a <- cFA
		x <- cJA
		z <- cJF
		f <- fA
		v <- vuln
		h <- hide
		s <- surv
		d <- DF
				
		if(stateName=="F0"){
			# dVal <- d*(Fo-F0) - a*F0* ( (s*f)/(x*(qE+1-s)) - s/x - (z*v*F0)/(x*(h+v+z*F0)) )
			dVal_grow <- d*(Fo-F0)
			dVal_cons <- - a*F0* ( (s*f)/(x*(qE+1-s)) - s/x - (z*v*F0)/(x*(h+v+z*F0)) )
			
		}
		if(stateName=="J0"){
			Q <- 1/(qE+1-s)
			# dVal <- f*s*J0*Q - x*J0^2*s*Q - s*J0 - (z*v*J0)/((h+v)/((d*Fo)/(d+a*s*J0*Q)) + z)
			dVal_grow <- f*s*J0*Q
			dVal_cons <- - x*J0^2*s*Q - s*J0 - (z*v*J0)/((h+v)/((d*Fo)/(d+a*s*J0*Q)) + z)
		}
		
		if(parts){
			dVal <- c(growth=unname(dVal_grow), consumption=unname(dVal_cons))
		}else{
			dVal <- unname(dVal_grow + dVal_cons)
		}

		dVal
		
	})
	
	return(out)
	
}


#' Change in F temporal dynamics wrt F
#' 
#' Define a function u=dF/dt, this gives du/dF. The function u is a function of F and parameters.
#' 
#' @param F abundance of planktivorous fish
#' @param pars named vector of parameters found in \code{\link{ecoStep}}. At a minimum, qE needs to be supplied. Optionally, other parameter arguments can be supplied here, but if they aren't, they are taken from the formals (defaults) of \code{\link{ecoStep}}.
#' 
#' @return numeric value that is du/dt
#' @examples
#' dU_dF(3.702216, pars=c(qE=0.5))
#' @export
dU_dF <- function(F, pars){
	
	parsF <- unlist(formals(ecoStep)[c("fA", "cJA", "cJF", "cFA", "vuln", "hide", "surv", "Fo", "DF")])
	
	if(missing(pars)){ # if a function requires qE, pars needs to be supplied (e.g., pars=c(qE=0.1))
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	
	out <- with(as.list(pars), {
		qE <- pars['qE']
		a <- cFA
		x <- cJA
		z <- cJF
		f <- fA
		v <- vuln
		h <- hide
		s <- surv
		d <- DF
	
		# d + (a*s*f)/(x*(qE+1-s)) + (2*z*v*a*F*(x*(h+v+z*F))-x*z^2*v*a*F^2)/(x*(qE+1-s))^2
		-d - (a*s*f)/(x*(qE+1-s)) + a*s/x + (2*a*z*v*F)/(x*(h+v+z*F)) + (x*a*v*z^2*F^2)/(x*(h+v+z*F))^2
		
	})
	
	return(unname(out))
	
}


#' Find Critical Values
#' 
#' Find the critical values of harvest in the trophic cascade model
#' 
#' @param parRange range of parameter (qE) values over which to search
#' @param stateRange range of state variable (J0) values over which to search
#' @param par_nGrid size of grid of parameter (qE) values
#' @param state_nGrid size of grid of state variable (J0) values
#' 
#' @return a vector of critical values of qE
#' @examples
#'  findCritHarv(par_nGrid=1E2, state_nGrid=1E3) # fast but very inaccurate
#' findCritHarv(par_nGrid=1E2, state_nGrid=5E4) # closer, but still off
#' @export
findCritHarv <- function(parRange=c(0.05,1.5), stateRange=c(0,2E3), par_nGrid=3E2, state_nGrid=5E5){
	stopifnot(par_nGrid>=30)
	
	# qE values to initially screen
	qE_vals <- seq(parRange[1], parRange[2], length.out=par_nGrid)
	
	# function to return roots at each qE (basically like 'getRoot' function)
	getStates <- function(qE, xR){
		rootSolve::uniroot.all(f=dFJ_dt_1state, interval=c(xR[1], xR[2]), pars=c(qE=qE), n=state_nGrid)
	}
	qE_states <- lapply(qE_vals, getStates, xR=stateRange)
	qE_nStates <- sapply(qE_states, length) # n states per qE
	
	# function to return qE value when number of states changes
	dStates <- function(x, qe){
		qe[(1+which(diff(x)!=0))]
	}
	critGuess0 <- dStates(qE_nStates, qe=qE_vals) # initial guess at qE critical values
	
	return(critGuess0)
	
	
	# nStates <- function(qE, xR){
# 		length(getStates(qE,xR))
# 	}
# 	critVal <- c()
# 	xR_new <- range(qE_states)
# 	crit_jitter <- diff(parRange)*0.1*c(-1,1)
# 	# crit_jitter <- c(-2, 2)*mean(diff(qE_vals)) # could probably change to c(-1,1)
# 	for(j in 1:length(critGuess0)){
# 		# xR_new <- range(qE_states[qE_vals==critGuess0[j]])
# 		qE_range_new <- crit_jitter+critGuess0[j]
# 		qE_vals_new <- seq(qE_range_new[1], qE_range_new[2], length.out=par_nGrid)
# 		qE_states_new <- sapply(qE_vals_new, nStates, xR=xR_new)
# 		critVal[j] <- dStates(qE_states_new, qe=qE_vals_new)
# 	}
#
# 	return(critVal)
	
}


#' Stability Classification
#' 
#' Find equilibria for a model of a trophic cascade
#' 
#' @param qE numeric scalar representing harvest rate; reasonable values might be between 0 and 2
#' @param pars a named numeric vector of optional additional parameter values to be passed to \code{func}
#' @param func function whose stability is to be classified; nominally this is \code{fs2D}, a wrapper for \code{\link{fishStep2D}}; currently this function only works with the default b/c \code{\link{getRoot_df}} hasn't been set up to accept other functions as an argument
#' @param Jvals,Mvals numeric vector of J0 (juvenile bass) or F0 (planktivorous fish) values over which to search for equilibria
#' @return a data.frame with a row for each equilibrium point found, and columns for the J0,F0 coordinates of those points, a character describing its stability classification, the trace (tr) and determinant (Delta) of the Jacobian matrix at that point, the 'discriminant' value (tr^2 - 4*delta), and the parameter values supplied through qE and pars.
#' @examples
#' stabClass(qE=0.65)
#' qEvec <- seq(0.15, 1.5, length.out=10)
#' lout <- lapply(qEvec, stabClass)
#' do.call(rbind, lout)
#' @export
stabClass <- function(qE, pars, func=fs2D, Jvals=seq(0,2E3,length.out=2), Fvals=seq(10,20, length.out=2)){
	requireNamespace("phaseR", quietly=TRUE)
	if(missing(pars)){pars <- NULL}
	if("qE"%in%names(pars)){pars <- pars[names(pars)!="qE"]}
	gridVals <- cbind(expand.grid(J0=Jvals, F0=Fvals), qE=qE)
	rs <- getRoot2D_df(gridVals, pars=pars)
	rs <- rs[complete.cases(rs),,drop=FALSE]
	urs <- rs[!duplicated(paste0(round(rs[,"J0"],4),round(rs[,"F0"],4))), c("qE","J0","F0"), drop=FALSE]

	stab <- function(x){
		st <- phaseR::stability(fs2D, y.star=x[-1,drop=FALSE], parameters=c(qE=qE,pars), summary=FALSE)
		# if(is.null(names(st$parameters))){names(st$parameters) <- "I"}
		o <- cbind(data.frame(J0=st$y.star[1], F0=st$y.star[2], classification=st$classification, Delta=st$Delta, discriminant=st$discriminant, tr=st$tr),as.list(st$parameters))
		rownames(o) <- NULL
		o$classification <- as.character(o$classification)
		return(o)
	}
	do.call('rbind', apply(urs, 1, stab))
}


#' Equations of Motion
#' 
#' Provides the motion of each state variable (change in that variable as a function of itself and parameters)
#' 
#' @param state0 numeric, value of a state variable (can be scalar or vector, not named)
#' @param stateVar character, name of state variable whose motion you want
#' @param pars named vector of parameters, needs to at least include qE
#' 
#' @return numeric scalar or vector indicating the rate of change of the state variable
#' 
#' @note These functions should essentially replace \code{\link{dFJ_dt_1state}}. The A motion function was not an option in \code{\link{dFJ_dt_1state}}, though. Both A motion and J motion are done without much algebra or substitution; however, getting the equation/ function for F motion requires doing some alebraic substitution to get the A state variable to cancel out. Dividing by A and letting it cancel probably involves some extra assumption that might make the F motion equation/ function invalid in certain areas of state-space.
#' 
#' @export
stateMotion <- function(state0, stateVar=c("A0","F0","J0"), pars){
	stateVar <- match.arg(stateVar)
	parsF <- unlist(formals(ecoStep)[c("fA", "cJA", "cJF", "cFA", "vuln", "hide", "surv", "Fo", "DF")])
	if(missing(pars)){ # if a function requires qE, pars needs to be supplied (e.g., pars=c(qE=0.1))
		pars <- parsF
	}else{
		pars <- c(pars, parsF[!names(parsF)%in%names(pars)])
	}
	
	out <- with(as.list(pars), {
		# Functions that calculate the state variable value
		Afun <- function(Jhat){ # gives A as a function of J and parameters
			surv*Jhat/(qE+1-surv)
		}
		Afun2 <- function(Fhat){ # gives A as a function of F and parameters, required some algebra first, though
			surv*fA/(cJA*(qE+1-surv)) - (cJF*vuln*Fhat)/(cJA*(hide+vuln+cJF*Fhat))
		}
		Ffun <- function(Ahat){ # gives F as a function of A and parameters
			DF*Fo/(cFA*Ahat+DF)
		}
		Jfun <- function(Ahat, Fhat){ # gives J as a function of A, F, and parameters
			fA*Ahat/(cJA*Ahat + surv + (cJF*vuln*Fhat)/(hide+vuln+cJF*Fhat))
			# f*Ahat / ( cJA * Ahat + cJF* vuln * Fhat / (hide + vuln + cJF * Fhat)   ) # what Carl has; dropping the surv term adds a lot more curvature to the A motion equation, so it's easier to produce alternative stable states. However, this was tested with the Carpenter et al. 2008 Ecology Letter parameterization, but that paper used the continuous + discrete version of the model.
		}
		
		# Functions that calculate the change in the state variable as a function of itself and parameters
		Amotion <- function(Ahat){ # gives dA/dt as a function of A and parameters
			Jhat <- Jfun(Ahat=Ahat, Fhat=Ffun(Ahat=Ahat))
			(surv)*Jhat - qE*Ahat - ((1-surv))*Ahat
			# (surv)*Jhat - qE*Ahat - (mA<-0)*Ahat # what carl has
		}
		Fmotion <- function(Fhat){ # gives dF/dt as a function of F and parameters
			Ahat <- Afun2(Fhat=Fhat) # note that Afun2 required doing some algebra first, and might be making unspecified assumptions to be valid; careful!
			DF*(Fo-Fhat) - cFA*Fhat*Ahat
		}
		Jmotion <- function(Jhat){ # gives dJ/dt as a function of J and parameters
			Ahat <- Afun(Jhat=Jhat)
			Fhat <- Ffun(Ahat=Ahat)
			Jpredloss <- (-cJA*Jhat*Ahat)-(cJF*vuln*Jhat*Fhat/(hide + vuln + cJF*Fhat) )
			(fA)*Ahat + Jpredloss - (surv)*Jhat
		}
		
		# select the appropriate motion function, and run it on the provided state value
		# motion <- switch(stateVar, "A0"=match.fun("Amotion"), "F0"=match.fun("Fmotion"), "J0"=match.fun("Jmotion"))
		motion <- switch(stateVar, "A0"=Amotion, "F0"=Fmotion, "J0"=Jmotion)
		motion(state0)
	})
	return(unname(out))
}


