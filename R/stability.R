
#' Claculate the Roots of F
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



#' Temporal dynamics of F
#' 
#' Temporal dynamics of F, but as a function of F and parameters (not other state variables)
#' 
#' @param F abundance of planktivorous fish
#' @param pars named vector of parameters found in \code{\link{ecoStep}}. At a minimum, qE needs to be supplied. Optionally, other parameter arguments can be supplied here, but if they aren't, they are taken from the formals (defaults) of \code{\link{ecoStep}}.
#' 
#' @return numeric value indicating temporal rate of change of F
#' @examples
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.5)) # find equilibria when bass are rare
#' dF_dt_1state(F0=17.890184, pars=c(qE=0.5)) # F0 set to equilibrium when bass start as rare
#' 
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.5)) # find equilibria when bass are abundant
#' dF_dt_1state(F0=0.06712064, pars=c(qE=0.5)) # F0 set to equilibrium when bass start as abundant
#' 
#' getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.65)) # find equilibria when bass are rare
#' dF_dt_1state(F0=100, pars=c(qE=0.65)) # F0 set to equilibrium when bass start as rare # WTF?!
#' fishStep(X=c(A0=6.275129e-18, F0=100, J0=1.443280e-17))
#' 
#' getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.65)) # find equilibria when bass are abundant
#' dF_dt_1state(F0=0.09136212, pars=c(qE=0.65)) # check planktivore equation for rate of change
#' fishStep(X=c(A0=364.51517642, F0=0.09136212, J0=838.38490576)) # re-examine rates of change of fish near equilibrium
#' @export
dF_dt_1state <- function(F0, pars){
	
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
		
		# d*(Fo-F) - a*F* ( (s*F)/(x*(qE+1-s)) - (z*v*F)/(x*(h+v+z*F)) ) # wrong
		# d*(Fo-F) - a*F* ( (s*f)/(x*(qE+1-s)) - s/x - (z*v*F)/(x*(h+v+z*F)) )
		d*(Fo-F0) - a*F0* ( (s*f)/(x*(qE+1-s)) - s/x - (z*v*F0)/(x*(h+v+z*F0)) )
		
		# Acalc = (s*J)/(qE+1-s)
		# Fcalc = (d*Fo)/(d + a*A)
		# Jcalc = (f*A)/(x*A + s + (z*v*F)/(h+v+z*F))
		# A_JinA = (s*f)/(x*(qE+1-s)) - s/x - (z*v*F)/(x*(h+v+z*F)) # value for A assuming all states at equilibrium
		# A_JinA2 = (s*Jcalc)/(qE+1-s)
		# A_JinA3 = (s*((f*A)/(x*A + s + (z*v*F)/(h+v+z*F))))/(qE+1-s)
		#
		#
		# A*(qE+1-s)*(x*A + s + (z*v*F)/(h+v+z*F)) == s*f*A
		# (qE+1-s)*(x*A + s + (z*v*F)/(h+v+z*F)) == s*f
		# x*A + s + (z*v*F)/(h+v+z*F) == (s*f)/(qE+1-s)
		# A == (s*f)/(x*(qE+1-s)) - s/x - (z*v*F)/(x*(h+v+z*F))
		#
		# d*(Fo-F) - a*f*A
		
	})
	
	return(unname(out))
	
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








