setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
#####################IVQR GMM CODE############################################################
##############################################################################################
##############################################################################################
#############################################################################################
# Smoothing function Itilde() and its derviative Itilde'(); Vector input supported
Itilde.KS17 <- function(u) {  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0)) }
Itilde.deriv.KS17 <- function(u) { ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0) }
#########################
gmmq <- function(tau, va, Xt, lX, Z, Lagphi, id, time, Lambda=function(va, Xt, lX, Lagphi, b) va-Xt%*%b[1:(length(b)-1)]-b[length(b)]*(Lagphi-lX%*%b[1:(length(b)-1)]), 
  Lambda.derivative=function(va, Xt, lX, Lagphi, b) va-Xt%*%b[1:(length(b)-1)]-b[length(b)]*(Lagphi-lX%*%b[1:(length(b)-1)]), h=0, VERBOSE=FALSE, b.init=0) {
  #Sample Size
  dZ <- ncol(Z)
  n <- length(va)
  # Set smoothing function Itilde() and its derviative Itilde'()
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
  #Objective Function
  obj.fn <- function(b,h) { #objective function
    L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b), ncol=1)
    return(colMeans(Z*array(data=Itilde(-L/h)-tau, dim=dim(Z))))
  }
  if (!is.null(Lambda.derivative)) {
    jac.fn <- function(b,h) {
      L <- Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b)
      Lp <- Lambda.derivative(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b)
      return((t(Z)%*%(Lp*array(data=Itilde.deriv(-L/h),dim=dim(Z))*(-1/h)))/n)
    }
  }
  # Compute estimator: if h=0, find smallest h newtonsys can use; else use user's h (if big enough)
  if (h==0) {
    h.init <- h.cur <- 0.001 #was: .Machine$double.eps * 1e3
    hfac <- 100
  } else {
    h.init <- h.cur <- h
    hfac <- 10
  }
  h.lo <- h
  h.hi <- NA
  b.best <- matrix(data=b.init,nrow=dZ,ncol=1)
  MAXITER <- 400;  H.MAX <- 1e10 
  # Loop: h<=h.lo<h.cur<h.hi, where h.hi can be solved by newtonsys but h.lo cannot
  if (VERBOSE) cat(sprintf('n=%d, tau=%g, dZ=%d, desired h=%g\n',n,tau,dZ,h))
  last.error <- NULL
  while (h.cur<=H.MAX && h.cur>.Machine$double.eps*1e2 && 
         (is.na(h.hi) || (h.hi/h.lo)>1.4)) {
    if (VERBOSE) cat(sprintf("h.cur=%g: ",h.cur))
    if (is.null(Lambda.derivative)) J <- NULL else J <- function(b)jac.fn(b,h.cur)
  	soln <- tryCatch(newtonsys(Ffun=function(b) obj.fn(b,h.cur), Jfun=J, x0=b.best, maxiter=MAXITER),
  				  warning=function(w) NA,
  				  error=function(e) {return(NA)} #warning("(above) error from running pracma::newtonsys."); stop(e) } #NA
  				)
  	exitOK <- !(is.na(soln[1]) || soln$niter==MAXITER)
  	if (!exitOK) {
  	  last.error <- soln[2]
  	  if (VERBOSE) {
  	    if (is.na(soln[1])) cat(sprintf("no solution, soln=NA\n")) else cat(sprintf("no solution, soln$niter=%d\n",soln$niter))
  	  }
  	  h.lo <- h.cur
  	  if (is.na(h.hi)) { #keep going up till find h.cur that actually works
  	    h.cur <- h.cur*hfac
  	  } else { #try arithmetic [not geometric] mean of h.lo and h.hi
  	    h.cur <- (h.lo+h.hi)/2 #sqrt(h.lo*h.hi)
  	  }
  	} else { #newtonsys found a solution
  	  if (VERBOSE) cat(sprintf("Solution found!\n"))
  	  h.hi <- h.cur; b.best <- soln$zero; h.lo <- h
  	  if (h.cur==h) { #found solution w/ user-requested h
  	    break
  	  } else if (h.hi/h.init>2.5 || h==0) { 
  	    h.cur <- (h.lo+h.hi)*(2/3)
  	  } else h.cur <- h.init
  	} 
  }
  if (h.cur>H.MAX) {
    warning("(above) error from pracma::newtonsys")
  }
  return(list(b=b.best,h=h.hi))
}































