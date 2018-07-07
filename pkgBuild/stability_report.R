#' ---
#' title: "Stability analysis of a trophic cascade model"
#' author: "Ryan Batt"
#' date: "2018-06-16"
#' abstract: |
#'        Makings of a stability analysis of a trophic cascade model. Starting by focusing on fish (adult bass, juvenile bass, and a planktivore), though the model is easily extensible to include phytoplankton and zooplankton for a 5D model.
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     fig_caption: true
#'     theme: "readable"
#'     template: default
#'   pdf_document:
#'     toc: true
#'     toc_depth: 4
#'     template: latex-ryan.template
#'     fig_caption: yes
#' geometry: margin=1.0in
#' lineno: true
#' lineSpacing: false
#' titlesec: true
#' documentclass: article
#' placeins: true
#' ---


#+ report-setup, include=FALSE, echo=FALSE, cache=FALSE
# ==========================================
# = Record Time to Get Elapsed Time at End =
# ==========================================
t1 <- Sys.time()


# =================
# = Load Packages =
# =================
library(viridis)
library(phaseR)
library(rootSolve)
library(bs.tipping)

# Report
library(knitr)
library(rmarkdown)


# ================
# = Report Setup =
# ================
doc_type <- c("html", "pdf")[1]
table_type <- c("html"="html", "pdf"="latex")[doc_type]
options("digits"=3) # rounding output to 4 in kable() (non-regression tables)
o_f <- paste(doc_type, "document", sep="_")

# render!
# rmarkdown::render(
# 	"~/Documents/School&Work/epaPost/bs-tipping/pkgBuild/stability_report.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/bs-tipping/pkgBuild/',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'stability_report/', 
	cache.path='stability_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=FALSE,
	autodep=TRUE,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


#' #Setup
#+ setup
qE <- 0.65
qE_end <- 0.1

dt <- 0.01
nYears <- 1000

noise_coeff <- c(0.01, 0.01, 0.01)


#' #Simulation
#+ simulation

# set initial values to equilibrium (roots)
X_init <- pmax(bs.tipping::getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=qE)), 1E-2)

# simulate an example of fish
qE_vec <- seq(qE, qE_end, length.out=nYears/dt)
stateMat <- matrix(NA, nrow=nYears/dt, ncol=3, dimnames=list(NULL, c("A0","F0", "J0")))
stateMat[1,] <- unlist(X_init)

for(j in 2:nrow(stateMat)){ # iterate through time steps
	state <- stateMat[j-1,]
	dState_dt <- fishStep(X=state, pars=c(qE=(qE_vec[j])))

	# Euler Method Approximation
	dState <- dState_dt*dt + rnorm(3, sd=c(noise_coeff[1], noise_coeff[2], noise_coeff[3]))*dt 
	eulerState <- pmax(state + dState, 1E-2)
	stateMat[j,] <- eulerState # euler
}

stateMat <- cbind(time=seq(0, nYears-dt, by=dt), qE=qE_vec, stateMat)

#' ##Plot Simulation
#+ plot-simulation, fig.width=3.5, fig.height=6, fig.cap="Simulated time series of adult bass (A0), planktivores (F0), and juvenile bass (J0) over a gradient of harvesting of adult bass (qE)."
par(mfrow=c(3,1), mar=c(2, 2, 0.75, 0.25), ps=8, cex=1, mgp=c(1, 0.25, 0), tcl=-0.15)
plot(stateMat[,c("time","A0")], type='l')
plot(stateMat[,c("time","F0")], type='l')
plot(stateMat[,c("time","J0")], type='l')
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Rearrange equations to represent stability in 1D
#' Start with the equations for the fish dynamics:  
#' \begin{array}{llr}
#' \dot{A} &= sJ - qEA - (1-s)A &(1) \\
#' \dot{F} &= d(F_o - F) - c_{FA}FA &(2)\\
#' \dot{J} &= fA - c_{JA}JA - \frac{c_{JF}vJF}{h+v+c_{JF}F} - sJ &(3)
#' \end{array}
#' 
#' 
#' Then, set $\dot{A}$ and $\dot{J}$ to 0, and solve for $A$ and $J$, respectively:  
#' \begin{array}{llr}
#' A&= \frac{sJ}{qE+1-s} &(4)\\[10pt]
#' J&=\frac{fA}{xA+s+\frac{zvF}{h+v+zF}} &(5)\\[10pt]
#' \end{array}
#' 
#' Substitute Eq5 into Eq4, and solve for $A$:  
#' \begin{array}{llr}
#' A&=\frac{sf}{x(qE+1-s)} - \frac{s}{x} - \frac{zvF}{x(h+v+zF)} &(6)
#' \end{array}
#' 
#' Substitute Eq6 into Eq2, giving us the dynamics of $F$ as a function of $F$ and parameters:  
#' \begin{array}{llr}
#' \dot{F} &= d(F_o-F) - aF(\frac{sf}{x(qE+1-s)} - \frac{s}{x} - \frac{zvF}{x(h+v+zF)}) &(7)
#' \end{array}
#'   
#' Eq7 is what is used in `dF_dt_1state`.  
#' 
#' **Alternatively**, we can rearrange the equations to perform a different seat of substituions to solve for $\dot{J}$ as a function of $J$ and parameters (juvenile bass):  
#' \begin{array}{llr}
#' Q &\equiv 1/(qE + 1 - s) \\
#'  \dot{J} &= fsJQ - xsQJ^2 - sJ - \frac{zvJ}{\frac{(h+v)}{(dF_o)(d+asQJ)^{-1}} + z} &(8)
#' \end{array}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Plots of dState/dt vs State
#' ##dF/dt vs F
#+ plot-dFdt-vs-F, fig.width=3.5, fig.height=3.5, fig.cap="dF/dt vs F. Where the line intersects 0, indicates a value of planktivore abundance that is an equilibrium. Is contingent upon the value of parameters, such as qE. The value of qE used here is the value used for the first time step of the simulation."
# check rate equation ... should be 0 at root
F_grad <- seq(-2, 110, length.out=100)
par(mar=c(2,2,0.75,0.5), mgp=c(1,0.25,0), tcl=-0.15, ps=8, cex=1)
plot(F_grad, dFJ_dt_1state(State0=F_grad, pars=c(qE=qE), stateName="F0"), type='l', ylab="dF/dt", xlab='F')
mtext(paste0('qE = ', round(qE,2)), line=-0.1, adj=0.05, font=2)
abline(h=0, lty=2)
#' This plot of the potential looks wrong to me. At a harvest rate of qE=0.65, F should have an equilibrium near F=100, but this plot shows a quadratic where the change in F is getting more and more positive past ~10. Furthermore, 10 looks like a saddle, not a stable node.  
#' 
#' ##dJ/dt vs J
#+ plot-dJdt-vs-J, fig.width=3.5, fig.height=3.5, fig.cap="dJ/dt vs J Where the line intersects 0, indicates a value of juvenile bass abundance that is an equilibrium. Is contingent upon the value of parameters, such as qE. The value of qE used here is the value used for the first time step of the simulation."
# check rate equation ... should be 0 at root
J_grad <- seq(0, 850, length.out=1E3)
par(mar=c(2,2,0.75,0.5), mgp=c(1,0.25,0), tcl=-0.15, ps=8, cex=1)
plot(J_grad, dFJ_dt_1state(State0=J_grad, pars=c(qE=qE), stateName="J0"), type='l', ylab="dJ/dt", xlab='J')
mtext(paste0('qE = ', round(qE,2)), line=-0.1, adj=0.05, font=2)
abline(h=0, lty=2)

#' ##dJ/dt vs J for various qE
#+ plot-dJdt-vs-J-qE, fig.width=3.5, fig.height=3.5, fig.cap="dJ/dt vs J Where the line intersects 0, indicates a value of juvenile bass abundance that is an equilibrium. Is contingent upon the value of parameters, such as qE. Thus, the dJ/dt vs J curve is plotted for several values of qE, providing a visualization of how the stability landscape changes with harvest rate."
# check rate equation ... should be 0 at root
J_grad <- seq(-1, 850, length.out=1E3)
# qE_vals <- seq(qE, qE_end, length.out=4)
qE_vals <- seq(1.25, 0.25, length.out=8)
dFJ_dt_1state_Jwrap <- function(X){dFJ_dt_1state(State0=J_grad, pars=c(qE=X), stateName="J0")}
dJ_dt_qE <- lapply(qE_vals, FUN=dFJ_dt_1state_Jwrap)
names(dJ_dt_qE) <- paste0("qE",round(qE_vals,2))

par(mar=c(2,2,0.75,0.5), mgp=c(1,0.25,0), tcl=-0.15, ps=8, cex=1)
ylim <- range(dJ_dt_qE)
cols <- viridis(n=length(qE_vals))
plot(J_grad, dJ_dt_qE[[1]], type='l', ylab="dJ/dt", xlab='J', ylim=ylim, col=cols[1])
for(j in 2:length(qE_vals)){
	lines(J_grad, dJ_dt_qE[[j]], col=cols[j])
}
legend('bottomleft', lty=1, col=cols, legend=paste0("qE = ",round(qE_vals,2)), ncol=2, y.intersp=0.6, x.intersp=0.5)
abline(h=0, lty=2)
#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Sanity Check on dF/dt
#' ##Part 1
#' The code below is taken from the examples in the help file for `dF_dt_1state`.  
#'   
#+ sanity-check1, results='markup'
getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.5)) # find equilibria when bass are abundant
dFJ_dt_1state(State0=0.06712064, pars=c(qE=0.5), stateName="F0") # F0 set to equilibrium when bass start as abundant
#' Check -- results make sense, the `getRoot` function finds that when bass start off super abundant, planktivores should stabilize near `0`. Using the `dF_dt_1state` function, we find that the change in planktivore abundance per unit time is very near `0` when planktivore abundance is near `0`. Great.  
#'   
#' ##Part 2
#+ sanity-check2, results='markup'
getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.5)) # find equilibria when bass are rare
dFJ_dt_1state(State0=17.890184, pars=c(qE=0.5), stateName='F0') # F0 set to equilibrium when bass start as rare
#' Check -- again, we find that according to `dF_dt_1state`, $dF/dt$ is near `0` when fish abundance is set near to equilibrium. This time the equilibrium value for planktivores was higher because the bass starting point was low. Right. Great.  
#'   
#' ##Part 3
#+ sanity-check3, results='markup'
getRoot(c(A0=1000, F0=1, J0=1000), pars=c(qE=0.65)) # find equilibria when bass are abundant
dFJ_dt_1state(State0=0.09136212, pars=c(qE=0.65), stateName="F0") # check planktivore equation for rate of change
fishStep(X=c(A0=364.51517642, F0=0.09136212, J0=838.38490576)) # re-examine rates of change of fish near equilibrium
#' Here we are repeating **Part 1**, but we've increased `qE` from `0.5` to `0.65`. I'm also cross-checking the dF/dt values from `dF_dt_1state` with those reported by `fishStep`, which is what `getRoot` uses to find equilibria. With the slightly higher harvest and high starting values for bass, we again find that planktivores should stabilize near `0`, and both `dF_dt_1state` and `fishStep` indicate that dF/dt is very small when setting F to this near-`0` equilibrium value. It's a little annoying that these two functions don't report the *exact* same value for dF/dt, but in general, this seems to check out. So far, things make sense.  
#'   
#' ##Part 4
#+ sanity-check4, results='markup'
getRoot(c(A0=1, F0=1, J0=1), pars=c(qE=0.65)) # find equilibria when bass are rare
dFJ_dt_1state(State0=100, pars=c(qE=0.65), stateName="F0") # F0 set to equilibrium when bass start as rare # WTF?!
fishStep(X=c(A0=6.275129e-18, F0=100, J0=1.443280e-17)) # but this one says that dF/dt should be ~0!!
#' Here we again use the `qE=0.65` harvest rate, with low starting values for bass. But instead of `dF_dt_1state` reporting a near-`0` value for the rate of change in planktivores, it's reporting a very large positive value. This last example is what is making me scratch my head.   
#'   
#' Update: apparently **Part 4** (above) doesn't work properly because in the algebra we dvide by A in a place that implies the assumption that A≠0. So when we try to use the equation when A=0, it doesn't work well.  
#' 
#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Stability Analysis
#' ##Critical Values
#+ critVals
(critVals <- findCritHarv()) # 0.24 and 1.22
#' These critical values were found numerically using the 1D model.  
#'   
#' ##Stability Classification
#+ stabClass, results="markup"
qEvec <- c(0, critVals[1]-0.2, critVals[1], critVals[2]-0.1, critVals[2]-0.05, critVals[2])
lout <- lapply(qEvec, stabClass)
do.call(rbind, lout)
#' I'm a little confused, because there is no place I could find with 2 stable nodes and 1 saddle point.  
#'   
#' ##Phase Portrait
#+ figure-phasePortrait, fig.width=6, fig.height=6, fig.cap="**Figure 1.** A phase portrait of the system for varying values of harvest (qE). The vector field (indicating the direction and speed that the system moves through phase space at that point) is represented by gray arrows. Nullclines are represented red and blue lines, indicating where dJ/dt and dF/dt are equal to zero, respectively.  Trajectories starting at arbitrary initial points (open diamonds) and continuing the along the accompanying solid black line indicate how the system moves from the initial point through phase space for 20 years. Equilibria are indicated by points: solid filled circle is a stable node, an 'X' is a saddle point. An equilibrium occurs whereever the nullclines cross. The different panels correspond to different values of harvest (qE). "
qEvals <- rev(c(critVals[1]-0.2, critVals[1], critVals[2]-0.1, critVals[2]))
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=9, cex.axis=0.85)
for(j in 1:4){
	(bs.tipping::phasePortrait(qE=qEvals[j], pars=NULL, nFlow=10, nNull=300, t.end=20, addLeg=TRUE))
	mtext(paste0("qE = ",round(qEvals[j],2)), side=3, line=0, adj=0, font=2)
}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' 
#' ##Growth and Consumption Curves: dJ/dt vs J
#+ figure-grow-cons-J, fig.width=6, fig.height=6, fig.cap="**Figure.** Growth and consumption curves for juvenile bass (J). Different panels show the curves for varying values of harvest rate on adult bass (qE)."
# dFJ_dt_1state <- function(State0, pars, stateName=c("J0","F0"), parts=FALSE){
plot_growCons <- function(stateRange=c(0,999), pars=c(qE=1.0), stateName=c("J0","F0"), nGrid=100){
	stateName <- match.arg(stateName)
	stopifnot(length(stateRange)==2)
	
	names(stateRange) <- c("from", "to")
	grid_args <- c(as.list(stateRange), length.out=nGrid)
	stateGrid <- do.call(seq, grid_args)
	
	rates <- t(sapply(stateGrid, dFJ_dt_1state, pars=pars, stateName=stateName, parts=TRUE))
	state_rates <- data.table(state=stateGrid, rates)
	
	ylim <- state_rates[,range(c(growth, consumption))] #range(state_rates[,c("growth", "consumption")])
	state_rates[,plot(state,growth, col='blue', type='l', ylim=ylim, xlab="", ylab="")]
	state_rates[,lines(state,consumption, col='red')]
	
}
qEvals <- rev(c(critVals[1]-0.2, critVals[1], critVals[2]-0.1, critVals[2]))
par(mfrow=c(2,2), mar=c(2,2,1,0.5), mgp=c(1,0.25,0), tcl=-0.15, cex=1, ps=9, cex.axis=0.85)
for(j in 1:4){
	
	plot_growCons(pars=c(qE=qEvals[j]))
	
	mtext(paste0("qE = ",round(qEvals[j],2)), side=3, line=0, adj=0, font=2)
	mtext("dJ/dt", side=2, line=0.85)
	mtext("J", side=1, line=0.85)
	legend("topleft", legend=c("growth", "consumption"), col=c("blue","red"), lty=1)
}


#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Session Info
#+ sessionInfo, results='markup'
difftime(Sys.time(), t1) # how long it took to run these models/ produce this report
Sys.time()
sessionInfo()
