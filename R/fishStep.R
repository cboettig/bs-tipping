# #' GAMMA
# #'
# #' Phytoplankton growth as a function of depth, photosynthesis-irradiance parameters, and phytoplankton abundance
# GAMMA <- function(z,Pbar) {
# 	Iz <- I0*exp(-z*(eps0+epsP*Pbar))
# 	rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
# }

#' Fish Step
#' 
#' Compute the rate of change of the three fish state variables given initial conditions and parameters
#' 
#' @param X vector of length 3 with names A0, F0, and J0.  Indicates starting value of state variables
#' @param pars a length 1 vector with name qE, indicating harvesting rate of adults
#' @param fA  Fecundity of adult piscivore (2 in OLD)
#' @param cJA Density dependent mortality rate of juveniles
#' @param cJF Consumption of juveniles by planktivores
#' @param cFA Consumption of planktivores by adult piscivores
#' @param vuln Vulnerability coefficient (this is "v" in eco lett table/ equations)
#' @param hide Hiding coefficient (this is "h" in eco lett table/ equations)
#' @param surv Overwinter survivorship of adults
#' @param Fo Refuge density of planktivores 
#' @param DF Diffusion parameter for planktivores
#' 
#' @return named vector of length 3, with names corresponding to juvenile bass (J0), adult bass (A0), and sunfish (F0) abundances.
#' @examples
#' fishStep(X=c(A0=1.5,F0=25,J0=2), pars=c(qE=0.1))
#' @export

# Phytoplankton Parameters
# I0 <- 300  # Surface irradiance, microEinsteins m-1 s-1
# # P-I curve parameters, Follows et al.
# k_sat <- 0.012 # per microEinsteins m-1 s-1
# k_inh <- 0.004 # per microEinsteins m-1 s-1 (nominal 0.004, range 0.001-0.007)
# # Light extinction parameters, Carpenter et al. L&O 1998
# DOC <- 5  # Assumed DOC conc, g m-3
# eps0 <- 0.0213 + 0.0514*DOC  # Baseline, per m (nominal DOC par 0.0514)
# epsP <- 0.0177   # Phytoplankton, per m2 (mg phosphorus)-1
#
# # Derived parameter from Follows et al.
# Fmax <- ((k_sat + k_inh)/k_sat)*exp(-(k_inh/k_sat)*log(k_inh/(k_sat+k_inh)))
#
# # Information for depth integration
# Zmix <- 4  # Mixed layer depth, m
# nZ <- 10  # Steps for vertical integration
# dZ <- Zmix/nZ
# Zvec <- seq(0,Zmix,by=dZ)
#
# rP <- 3  # Phytoplankton growth parameter per unit phosphorus load
# Load <- 0.6  # Daily phosphorus load
# mP <- 0.1  # Phytoplankton daily mortality

# Zooplankton Parameters
# Ho <- 1 #4    # Refuge biomass
# DH <- 0.5  # Diffusion parameter for herbivore (zooplankton)
# cHF <- 0.1  # Consumption rate by planktivore
# alf <- 0.3  # Conversion efficiency of consumed phytoplankton to zooplankton
# cPH <- 0.25  # Consumption rate of phytoplankton by zooplankton

# Fish Parameters
# fA <- 2  # Fecundity of adult piscivore (2 in OLD)
# cJA <- 1E-3 #0.1  # Density dependent mortality rate of juveniles
# cJF <- 0.5  # Consumption of juveniles by planktivores
# cFA <- 0.3  # Consumption of planktivores by adult piscivores
# vuln <- 1 #80  # Vulnerability coefficient (this is "v" in eco lett table/ equations)
# hide <- 8 #80  # Hiding coefficient (this is "h" in eco lett table/ equations)
# surv <- 0.5 #0.6  # Overwinter survivorship of adults
# Fo <- 100  # Refuge density of planktivores  # 100 in OLD
# DF <- 0.1 #0.09  # Diffusion parameter for planktivores
# sigma <- 0.1 #0.05  # SD of additive noise for planktivores (0.1 in May 07)
# A2biom <- 0.2  # Convert A to kg / ha
# J2biom <- 0.05  # Convert J to kg / ha
# F2biom <- 1  # Convert F to kg / ha

fishStep <- function(X, pars=c(qE=1), fA=2, cJA=1E-3, cJF=0.5, cFA=0.3, vuln=1, hide=8, surv=0.5, Fo=100, DF=0.1){
	# Ho=1, DH=0.5, cHF=0.1, alf=0.3, cPH=0.25, sigma=0.1,  A2biom=0.2, J2biom=0.05, F2biom=1
	with(as.list(c(X,pars)),{
		# Fish dynamics
		#Arate <- (surv/nint)*J0 - qE*A0 - ((1-surv)/nint)*A0
		Arate <- (surv)*J0 - qE*A0 - ((1-surv))*A0
		Frate <- DF*(Fo-F0) - cFA*F0*A0
		Jpredloss <- (-cJA*J0*A0)-(cJF*vuln*J0*F0/(hide + vuln + cJF*F0) ) # Note this is negative
	
		#Jrate <- (fA/nint)*A0 + Jpredloss - (surv/nint)*J0 
		Jrate <- (fA)*A0 + Jpredloss - (surv)*J0
		# A1 <- A0 + (Arate*dt)   # Update A numerically
	# 	# if(qE<0 & A1<1){A1 <- max(A1,7)} # RDB
	# 	F1 <- F0 + (Frate*dt) + (sigma*NoiseF*dtZ)
	# 	J1 <- J0 + (Jrate*dt)
	# 	A1 <- max(A1,0.1)  # Force A1 greater than 0.1
	# 	F1 <- max(F1,0.1)  # Force F greater than 0.1
	# 	J1 <- max(J1,0.1)  # Force J greater than 0.1
	
		# # Zooplankton dynamics
		# Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
		# H1 <- H0 + (Hrate*dt) + (sigma*NoiseH*dtZ)
		# H1 <- max(H1,0.1)  # Force H greater than 0.01
		#
		# # Phytoplankton dynamics
		# Pbar <- P0   # Set P value for vertical integration
		# gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
		# gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
		# Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0)
		# P1 <- P0 + (Prate*dt) + (sigma*NoiseP*dtZ)
		# P1 <- max(P1,0.1)  # Force P greater than 0.1

		# Construct output vector
		# simOut <- c("At"=A1, "Ft"=F1, "Jt"=J1, "Ht"=H1, "Pt"=P1)
		# SimList <- list(A1,F1,J1,H1,P1)
		dXdt <- c(dA_dt=Arate, dF_dt=Frate, dJ_dt=Jrate)
		return(dXdt)
		# pmax(dXdt, 0.1-c(A0,F0,J0)) # use this to make it so that the rate doesn't set a state variable below 0.1 ... but use in simulation, not in this function
		
	})
	
}

#' Two-dimensional fish time step
#' 
#' Rate of change in planktivores and juvenile bass per time step
#' 
#' @param X length 2 named vector of state variables. Names are should be J0 (juveniles) and F0 (planktivores)
#' @param pars a named vector of model parameters; in partiular, qE is harvest
#' @param ... Optional model parameters (e.g., see arguments in \code{\link{ecoStep}}), or other arguments to be passed to future methods. 
#' 
#' @return length-2 named vector (J0 and F0) indicating the rate of change per unit time in the abundances of juvenile bass (J0) and planktivores (F0).
#' 
#' @note \code{pars} and \code{...} are largely redundant. Their separation is merely for compatability with other functions, or for organizational purposes. The intended use is that parameters in \code{pars} correspond to parameters varied systematically to effect stability properties, whereas parameters specified through \code{...} are being altered for exploration purposes, but would generally not be supplied because the defaults in \code{\link{ecoStep}} are suitable.
#' 
#' @seealso \code{\link{fishStep}}, \code{\link{ecoStep}}
#' 
#' @examples
#' # set up parameter, initials, roots
#' p <- c(qE=0.65) # harvest rate
#' state1 <- c(A0=1, F0=1, J0=1) # initial values
#' state2 <- c(A0=1E3, F0=1, J0=1E3) # alternate initials
#' r1 <- getRoot(state1, pars=p)[c("J0","F0")] # first root
#' r2 <- getRoot(state2, pars=p)[c("J0","F0")] # second root (alt stbl states)
#' 
#' # test dState/dt at root
#' fishStep2D(r1, pars=p) # 0 at lower root
#' fishStep2D(r2, pars=p) # 0 at upper root
#' 
#' # test dState/dt near, but not at, root
#' fishStep2D(r1+rnorm(2), pars=p) # non-0 away from root
#' fishStep2D(r2+rnorm(2), pars=p) # non-0 away from root
#' 
#' # test dState/dt at r1 and r2, but different qE
#' fishStep2D(r1) # qE defaults to 1; r1 ~0 for qE=1 & =0.65
#' fishStep2D(r2) # r2 ~0 for qE=0.65, not for qE=1
#' @export
fishStep2D <- function(X, pars=c(qE=1), ...){
	parsF <- unlist(formals(ecoStep)[c("fA", "cJA", "cJF", "cFA", "vuln", "hide", "surv", "Fo", "DF")])
	
	dots <- list(...)
	
	parsFDots <- parsF
	
	if(length(dots)>0){
		dots_in_formals <- names(dots)%in%names(parsF)
		if(any(dots_in_formals)){
			pars_dots <- dots[dots_in_formals]
			parsFDots <- c(pars_dots, parsF[!names(parsF)%in%names(pars_dots)])
		}
	}
	
	if(missing(pars)){ # if a function requires qE, pars needs to be supplied (e.g., pars=c(qE=0.1))
		message("arg 'pars' missing, setting qE to 1")
		pars <- c(qE=1, parsFDots)
	}else{
		pars <- c(pars, parsFDots[!names(parsFDots)%in%names(pars)])
	}
	
	# for(j in 1:length(c(X,pars))){assign(names(c(X,pars))[j], c(X,pars)[j])}
	out <- with(as.list(c(X,pars)), {
		qE <- pars['qE']
		a <- cFA
		x <- cJA
		z <- cJF
		f <- fA
		v <- vuln
		h <- hide
		s <- surv
		d <- DF
		
		Q <- 1/(qE+1-s)
		dF_dt <- d*(Fo - F0) - a*F0*(s*J0*Q)
		dJ_dt <- f*(s*J0*Q) - x*J0*(s*J0*Q) - (z*v*J0*F0)/(h+v+z*F0) - s*J0
		
		dVal <- c(J0=unname(dJ_dt), F0=unname(dF_dt))
		dVal
		
	})
	# rm(list=c(c(X,pars), "qE", "a", "x", "z", "f", "v", "h", "s", "d", "Q", "dF_dt", "dJ_dt", "dVal"))
	
	# return(unname(out))
	return(out)
}


#' Wrapper for model when finding stability properties
#' 
#' Wrapper for 2D version of model for finding stability properties using \code{phaseR} package
#' 
#' @param t time, for compatibility with \code{phaseR} and \code{deSolve} and \code{rootSolve}
#' @param y state variables as in \code{\link{fishStep2D}}
#' @param parameters same as \code{pars} in \code{\link{fishStep2D}}; had to be called 'parameters' for compatibility with \code{phaseR}
#' @param ... other arguments to be passed to \code{\link{fishStep2D}}
#' @return length-1 list of length-2 named numeric vector indicating rate of change of state variables
#' @seealso \code{\link{fishStep2D}}
#' @export
fs2D <- function(t, y, parameters, ...){
	names(y) <- c("J0","F0")
	list(fishStep2D(y, pars=parameters, ...))
}
