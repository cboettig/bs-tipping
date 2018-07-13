#' ---
#' title: "Bifurcation or state tipping: assessing transition type in a model trophic cascade"
#' author: "Ryan Batt"
#' date: "2018-07-07"
#' abstract: |
#'        The primary purpose of this document is to sketch out the key equations, figures, and numerical results for a manuscript that focuses on assessing whether a model undergoes state-tipping or a bifurcation. Our approach emphasizes reducing the dimensionality of the original model to redefine manipulated state variables as parameters. Thus, when these redefined parameters are manipulated to induce a transition, because they parameters not states, the transition should be viewed as a bifurcation, not state tipping.
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
# 	"~/Documents/School&Work/epaPost/bs-tipping/pkgBuild/manuscript_report.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/epaPost/bs-tipping/pkgBuild/',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))
opts_chunk$set(
	fig.path = 'manuscript_report/', 
	cache.path='manuscript_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=FALSE,
	autodep=TRUE,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)



#' #System of Equations & Manipulation
#' Start with the equations for the ecosystem dynamics:  
#' \begin{array}{llr}
#' \dot{A} &= sJ - qEA - m_AA &(1) \\
#' \dot{F} &= D_F(F_o - F) - c_{FA}FA &(2)\\
#' \dot{J} &= fA - c_{JA}JA - \frac{c_{JF}vJF}{h+v+c_{JF}F} - sJ &(3)\\
#' \dot{H} &= D_H(H_R-H) + \alpha c_{HP}HP - c_{HF}HF &(4)\\
#' \dot{P} &= r_P L \gamma (I_0,P)P - mP - c_{PH}HP &(5)\\
#' \end{array}
#' 
#' 
#' Focusing on the fish, and setting each of the 3 fish differential equations to 0 and solving for the respective state variables:  
#' \begin{array}{llr}
#' A&= \frac{sJ}{qE+m_A} &(6)\\[10pt]
#' F&=\frac{D_FF_0}{D_F + c_{FA}A} &(7)\\[10pt]
#' J&=\frac{fA}{c_{JA}A+s+\frac{c_{JF}vF}{h+v+c_{JF}F}} &(8)\\[10pt]
#' \end{array}
#' 

#' 
#' **Begin Equations of Motion**  
#' Although our ultimate goal is to see if the equation of motion for $J$ or $F$ shows alternate states, for completeness we will also present the equation of motion for $A$, which is acquired by substituting the RHS of Eq7 for $F$ in Eq8, then substituting the resulting expression for $J$ in Eq1:  
#' \begin{array}{llr}
#' \dot{A} &=  sfA\Bigg(c_{JA}+s+c_{JF}v\bigg(\frac{(h+v)(D_F+c_{FA}A)}{D_F F_o}+c_{JF}\bigg)^{-1}\Bigg)^{-1} &(9)
#' \end{array}
#' 
#' To get the equation of motion for $F$, we begin by getting an expression for $A$ that depends only on $F$ and parameters, which is achieved by substituting Eq8 into Eq6 and solving for $A$:
#' \begin{array}{llr}
#' A&=\bigg({\frac{sfA}{c_{JA}A+s+\frac{c_{JF}vF}{h+v+c_{JF}F}}}\bigg)\bigg(qE+m_A\bigg)^{-1} &(10a)\\[15pt]
#' sfA &= A(qE+m_A)\big(c_{JA}A+s+\frac{c_{JF}vF}{h+v+c_{JF}F}\big) &(10b)\\[10pt]
#' sf &=(qE+m_A)\big(c_{JA}A+s+\frac{c_{JF}vF}{h+v+c_{JF}F}\big) &(10c)\\[10pt]
#' A&=\frac{sf}{c_{JA}(qE+m_A)} - \frac{s}{c_{JA}} - \frac{c_{JF}vF}{c_{JA}(h+v+c_{JF}F)} &(10d)
#' \end{array}
#' Going from Eq10b to Eq10c is achieved by dividing both sides of the equation by $A$, which implies that in Eq10c $A$ cannot be 0 when $\dot{A}$ is 0. We will return to the consequences of this assumption later as we assess our derived equation of motion for $F$. By subtituting the RHS of Eq10d for $A$ in Eq2, we arrive at the equation of motion for $F$:
#' \begin{array}{llr}
#' \dot{F} &= D_F(F_o-F) - c_{FA}F(\frac{sf}{c_{JA}(qE+m_A)} - \frac{s}{c_{JA}} - \frac{c_{JF}vF}{c_{JA}(h+v+c_{JF}F)}) &(11)
#' \end{array}
#'   
#' Eq10 is what is used in `dF_dt_1state`.  
#' 
#' The equation of motion for $J$ can be acquired by substituting the RHS of Eq7 for $F$ in Eq3 ,then exchanging $A$ in that resulting expression for the RHS of Eq6:  
#' \begin{array}{llr}
#' Q &\equiv 1/(qE + m_A) \\
#'  \dot{J} &= fsJQ - c_{JA}sQJ^2 - sJ - \frac{c_{JF}vJ}{\frac{(h+v)}{(D_FF_o)(D_F+c_{FA}sQJ)^{-1}} + c_{JF}} &(12)
#' \end{array}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   

#' #Setup
#+ setup
set.seed(42)
qE <- -0.001 #0.02
qE_end <- 0.001#/20 #-5#-0.025

dt <- 0.01 # 0.01
nYears <- 1020#/20

noise_coeff <- c(0.01, 0.01, 0.01, 0.01, 0.01)


#' #Simulation
#' ##Simulate Time series
#+ simulation
# set initial values to equilibrium (roots)
X_init <- pmax(bs.tipping::getRoot5D(c(A0=150, F0=50, J0=100, H0=100, P0=500), pars=c(qE=qE)), 1E-3)

# simulate an example of fish
qE_vec <- seq(qE, qE_end, length.out=nYears/dt)
stateMat <- matrix(NA, nrow=nYears/dt, ncol=5, dimnames=list(NULL, c("A0","F0","J0","H0","P0")))
stateMat[1,] <- unlist(X_init)

for(j in 2:nrow(stateMat)){ # iterate through time steps
	state <- stateMat[j-1,]
	dState_dt <- ecoStep(X=state, pars=c(qE=(qE_vec[j])))

	# Euler Method Approximation
	dState <- dState_dt*dt + rnorm(5, sd=c(noise_coeff))*dt
	eulerState <- pmax(state + dState, 1E-3)
	stateMat[j,] <- eulerState # euler
	# stateMat[j,] <- getRoot5D(stateMat[j-1,], pars=c(qE=qE_vec[j]))
}

tsMat <- cbind(time=seq(0, nYears-dt, by=dt), qE=qE_vec, stateMat)

#' ##Figure: 5-dimension simulation
#+ figure-5dsimulation, fig.width=3.5, fig.height=6, fig.cap="**Figure.** Simulated time series of adult bass (A0), planktivores (F0), and juvenile bass (J0) over a gradient of harvesting of adult bass (qE)."
par(mfrow=c(5,1), mar=c(1.5, 2, 0.75, 0.25), ps=8, cex=1, mgp=c(1, 0.25, 0), tcl=-0.15)
stateLabs <- c(
	A0="A (Adult Bass)",
	F0="F (Planktivores)",
	J0="J (Juvenile Bass)",
	H0="H ( Herbivores)",
	P0="P (Phytoplankton)"	
)
stateNames <- names(stateLabs)
for(i in 1:length(stateNames)){
	plot(tsMat[,"time"], tsMat[,stateNames[i]], type='l', ylab=stateLabs[stateNames[i]], xlab="")
}
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   


#' #Equations of Motion
#' ##Figure: A Motion
#+ figure-eqMotion-A, fig.width=3.5, fig.height=3.5, fig.cap="**Figure.** Equation of motion for A."
qevals <- c(0.25, 0.025, 0, -0.025)
par(mfrow=c(2,2), mar=c(1.85,1.85,0.75,0.25), mgp=c(1,0.15,0), tcl=-0.15, ps=9, cex=1)
for(i in 1:length(qevals)){
	curve(stateMotion(state0=x, pars=c(qE=qevals[i])), from=0, to=20, ylab="dA/dt", xlab="A")
	abline(h=0, lty=2)
	mtext(paste0("qE=",round(qevals[i],2)), side=3, line=0, adj=0.05, font=2)
}
#' Looks good, the second stable node appears when qE>0
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure: J Motion
#+ figure-eqMotion-J, fig.width=3.5, fig.height=3.5, fig.cap="**Figure.** Equation of motion for J."
par(mfrow=c(2,2), mar=c(1.85,1.85,0.75,0.25), mgp=c(1,0.15,0), tcl=-0.15, ps=9, cex=1)
for(i in 1:length(qevals)){
	curve(stateMotion(state0=x, stateVar="J0", pars=c(qE=qevals[i])), from=0, to=20, ylab="dJ/dt", xlab="J")
	abline(h=0, lty=2)
	mtext(paste0("qE=",round(qevals[i],2)), side=3, line=0, adj=0.05, font=2)
}
#' As with the adults, the second stable node appears when qE>0
#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure: F Motion
#+ figure-eqMotion-F, fig.width=3.5, fig.height=3.5, fig.cap="**Figure.** Equation of motion for F."
par(mfrow=c(2,2), mar=c(1.85,1.85,0.75,0.25), mgp=c(1,0.15,0), tcl=-0.15, ps=9, cex=1)
for(i in 1:length(qevals)){
	curve(stateMotion(state0=x, stateVar="F0", pars=c(qE=qevals[i])), from=0, to=20, ylab="dF/dt", xlab="F")
	abline(h=0, lty=2)
	mtext(paste0("qE=",round(qevals[i],2)), side=3, line=0, adj=0.05, font=2)
}
#' Unfortunately, this state variable only shows a trivial stable node and a saddle point. This is the state variable that requires algebra cancelling the 'A' state variable in order to get its equation of motion (without A canceling, you cannot get a dF/dt vs F type plot easily [as a function of F and parameters]). We think this result is invalid, as there must be some assumptions in the algebra that invalidate our desired interpretation.  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure: Combined A, J, F Motion
#+ figure-eqMotion-AJF, fig.width=5, fig.height=5, fig.cap="**Figure.** Motion of A, J, and F. Note that the equation of motion for F failed to exhibit alternate states. In formulating the equation for motion of F, there is an algebraic step that requires dropping from terms via division by A, which assume that A≠0. Therefore, states of F that correspond to A=0 are invalid. The equation of motion for A indicates that one of the stable nodes is at A=0, implying that the equation of motion for F is inaccurate in the corresponding region."
xlabs <- c("A","J","F") #ylabs <- c("dA/dt", "dJ/dt", "dF/dt")
ylabs <- paste0("d",xlabs,"/dt")
svars <- paste0(xlabs,"0")
par(mfrow=c(3,3), mar=c(1.85,1.85,0.75,0.25), mgp=c(1,0.15,0), tcl=-0.15, ps=9, cex=1)
for(s in 1:3){
	for(i in 2:length(qevals)){
		curve(stateMotion(state0=x, stateVar=svars[s], pars=c(qE=qevals[i])), from=0, to=20, ylab=ylabs[s], xlab=xlabs[s])
		abline(h=0, lty=2)
		mtext(paste0("qE=",round(qevals[i],2)), side=3, line=0, adj=0.05, font=2)
	}
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