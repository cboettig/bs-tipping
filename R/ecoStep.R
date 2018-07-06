# #' GAMMA
# #'
# #' Phytoplankton growth as a function of depth, photosynthesis-irradiance parameters, and phytoplankton abundance
# GAMMA <- function(z,Pbar) {
# 	Iz <- I0*exp(-z*(eps0+epsP*Pbar))
# 	rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
# }

#' Ecosystem Step
#' 
#' Compute the rate of change of the five ecosystem state variables given initial conditions and parameters
#' 
#' @param X vector of length 3 with names P0 (initial phytoplankton), H0 (initial herbivore [zooplankton]), A0 (adult piscivorous fish), F0 (planktivorous fish), and J0 (juvenile piscivores).  Indicates starting value of state variables
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
#' @param I0 surface irradiance, microEinsteins m-1 s-1
#' @param k_sat,k_inh photosynthesis-irradiance parameters
#' @param DOC,eps0,epsP Light extinction parameters from Carpenter et al. L&O 1998; DOC is dissolved organic carbon; eps are in units per m
#' @param Zmix depth of mixed layer in meters
#' @param nZ number of steps for vertical integration
#' @param rP phytoplankton growth parameter
#' @param Load daily phosphorus load
#' @param mP phytoplankton daily mortality
#' @param Ho refuge biomass
#' @param DH diffusion parameter for zooplankton (herbivore)
#' @param cHF consumption rate by planktivore (of zooplankton)
#' @param alf conversion efficiency of consumed phytoplankton to zooplankton
#' @param cPH consumption rate of phytoplankton by zooplankton
#' @param I0 surface irradiance
#' @param k_sat per microEinsteins m-1 s-1
#' @param k_inh per microEinsteins
#' @param DOC dissolved organic carbon, g m-3
#' @param eps0,epsP extinction coefficient, per m
#' 
#' @note \code{F} are referred to as planktivorous fish, but they do feed on the juvenile piscivores (\code{J}).
#' 
#' @return named vector of length 5, with names corresponding to juvenile bass (J0), adult bass (A0), and sunfish (F0) abundances.
#' @examples
#' ecoStep(X=c(A0=1.5,F0=25,J0=2, H0=1, P0=1), pars=c(qE=0.01))
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
# cJA <- 0.1 #1E-3   # Density dependent mortality rate of juveniles
# cJF <- 0.5  # Consumption of juveniles by planktivores
# cFA <- 0.3  # Consumption of planktivores by adult piscivores
# vuln <- 80 #1  # Vulnerability coefficient (this is "v" in eco lett table/ equations)
# hide <- 80 #8  # Hiding coefficient (this is "h" in eco lett table/ equations)
# surv <- 0.6 #0.5 #0.6  # Overwinter survivorship of adults
# Fo <- 200 #100  # Refuge density of planktivores  # 100 in OLD
# DF <- 0.09 # 0.1  # Diffusion parameter for planktivores
# sigma <- 0.1 #0.05  # SD of additive noise for planktivores (0.1 in May 07)
# A2biom <- 0.2  # Convert A to kg / ha
# J2biom <- 0.05  # Convert J to kg / ha
# F2biom <- 1  # Convert F to kg / ha

# ecoStep <- function(X, pars=c(qE=1), fA=2, cJA=1E-3, cJF=0.5, cFA=0.3, vuln=1, hide=8, surv=0.5, Fo=100, DF=0.1, Zmix=4, nZ=10, rP=3, Load=0.6, mP=0.1, Ho=1, DH=0.5, cHF=0.1, alf=0.3, cPH=0.25){
ecoStep <- function(X, pars=c(qE=0.01), fA=2, cJA=0.1, cJF=0.5, cFA=0.3, vuln=80, hide=80, surv=0.6, Fo=200, DF=0.09, Zmix=4, nZ=10, rP=3, Load=0.6, mP=0.1, Ho=4, DH=0.5, cHF=0.1, alf=0.3, cPH=0.25, I0=300, k_sat=0.012, k_inh=0.004, DOC=5, eps0=0.0213+0.0514*DOC, epsP=0.0177){
	with(as.list(c(X,pars)),{
		# Plankton Dynamics
		Fmax <- ((k_sat + k_inh)/k_sat)*exp(-(k_inh/k_sat)*log(k_inh/(k_sat+k_inh)))
		GAMMA <- function(z,Pbar) {
			Iz <- I0*exp(-z*(eps0+epsP*Pbar))
			rate <- (1/Fmax)*(1 - exp(-k_sat*Iz))*exp(-k_inh*Iz)
		}
		dZ <- Zmix/nZ
		Zvec <- seq(0,Zmix, by=dZ)
		
		
		# Fish dynamics
		#Arate <- (surv/nint)*J0 - qE*A0 - ((1-surv)/nint)*A0
		# Arate <- (surv)*J0 - qE*A0 - ((1-surv))*A0
		Arate <- (surv)*J0 - qE*A0 - (mA<-0.4)*A0
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
		Hrate <- DH*(Ho-H0) + alf*cPH*H0*P0 - cHF*H0*F0
		# H1 <- H0 + (Hrate*dt) + (sigma*NoiseH*dtZ)
		# H1 <- max(H1,0.1)  # Force H greater than 0.01
		#
		# # Phytoplankton dynamics
		Pbar <- P0   # Set P value for vertical integration
		gamvec <- GAMMA(Zvec,Pbar) # vertical sequence of light effect on growth
		gamI <- dZ*sum(gamvec)  # vertically integrated light effect on growth
		Prate <- (rP*Load*gamI*P0) - (mP*P0) - (cPH*H0*P0)
		# P1 <- P0 + (Prate*dt) + (sigma*NoiseP*dtZ)
		# P1 <- max(P1,0.1)  # Force P greater than 0.1

		# Construct output vector
		# simOut <- c("At"=A1, "Ft"=F1, "Jt"=J1, "Ht"=H1, "Pt"=P1)
		# SimList <- list(A1,F1,J1,H1,P1)
		dXdt <- c(dA_dt=Arate, dF_dt=Frate, dJ_dt=Jrate, dP_dt=Prate, dH_dt=Hrate)
		return(dXdt)
		# pmax(dXdt, 0.1-c(A0,F0,J0)) # use this to make it so that the rate doesn't set a state variable below 0.1 ... but use in simulation, not in this function
		
	})
	
}