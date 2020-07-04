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
gmmq <- function(tau, Y, mX, mlX, mZ, vphi, vlag.phi, Lambda=function(tau, theta, Y, mZ, mX, mlX, vphi, vlag.phi){Y-mX%*%theta}, 
  Lambda.derivative=function(tau, theta, h, Y, mZ, mX, mlX, vphi, vlag.phi){Y-mX%*%theta[1:(ncol(mX))]-theta[length(theta)]*(vlag.phi-mlX%*%theta[1:(ncol(mX))])}, h=0, VERBOSE=FALSE, theta.init=0) {
  #Sample Size
  dZ <- ncol(mZ)
  n <- dim(Y)[1]
  # Set smoothing function Itilde() and its derviative Itilde'()
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
  #Objective Function
  obj.fn <- function(theta,h) { #objective function
    L <- matrix(Lambda(tau=tau, theta=theta, Y=Y, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi), ncol=1)
    return(colMeans(mZ*array(data=Itilde(-L/h)-tau, dim=dim(mZ))))
  }
  if (!is.null(Lambda.derivative)) {
    jac.fn <- function(theta,h) {
      L <- Lambda(tau=tau, theta=theta, Y=Y, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi)
      Lp <- Lambda.derivative(tau=tau,theta=theta, Y=Y, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi)
      return((t(mZ)%*%(Lp*array(data=Itilde.deriv(-L/h),dim=dim(mZ))*(-1/h)))/n)
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
  theta.best <- matrix(data=theta.init,nrow=dZ,ncol=1)
  MAXITER <- 400;  H.MAX <- 1e10 
  # Loop: h<=h.lo<h.cur<h.hi, where h.hi can be solved by newtonsys but h.lo cannot
  if (VERBOSE) cat(sprintf('n=%d, tau=%g, dZ=%d, desired h=%g\n',n,tau,dZ,h))
  last.error <- NULL
  while (h.cur<=H.MAX && h.cur>.Machine$double.eps*1e2 && 
         (is.na(h.hi) || (h.hi/h.lo)>1.4)) {
    if (VERBOSE) cat(sprintf("h.cur=%g: ",h.cur))
    if (is.null(Lambda.derivative)) J <- NULL else J <- function(theta)jac.fn(theta,h.cur)
  	soln <- tryCatch(newtonsys(Ffun=function(theta) obj.fn(theta,h.cur), Jfun=J, x0=theta.best, maxiter=MAXITER),
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
  	  h.hi <- h.cur; theta.best <- soln$zero; h.lo <- h
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
  return(list(theta=theta.best,h=h.hi))
}
# Kernel functions from Andrews (1991) eqn (2.7)
QS.fn <- function(x) ifelse(x==0,1,(25/(12*pi^2*x^2))*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)))
Bartlett.fn <- function(x) ifelse(abs(x)>=1,0,1-abs(x))
uniform.fn <- function(x) ifelse(abs(x)>=1,0,1) #a.k.a. "Truncated" (Andrews 1991)
# Optimal bandwidths from (5.2) (or 6.2) in Andrews (1991)
QS.ST.fn <- function(alpha2,n) 1.3221*(alpha2*n)^(1/5)
Bartlett.ST.fn <- function(alpha1,n) 1.1447*(alpha1*n)^(1/3)
uniform.ST.fn <- function(alpha2,n) 0.6611*(alpha2*n)^(1/5) #footnote 5
# (6.4) in Andrews (1991)
alpha2.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/(1-rhos)^8) / sum(ws*sigmas^4/(1-rhos)^4)
alpha1.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/((1-rhos)^6*(1+rhos)^2)) / sum(ws*sigmas^4/(1-rhos)^4)
#
LRV.est.fn <- function(tau, Y, mX, mlX, mZ, vphi, vlag.phi, Lambda, theta.hat, Itilde, h, structure=c('iid','ts','cluster'), cluster.X.col, LRV.kernel=c('QS','Bartlett','uniform'), LRV.ST=NA, VERBOSE=FALSE) {
  # if (missing(structure) || !is.character(structure)) stop("Argument structure must be 'iid' or 'ts' or 'cluster'")
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  n <- dim(mZ)[1]
  if (structure %in% c('iid','ts')) {
    if (structure=='iid') {
      LRV.kernel <- 'uniform'; LRV.lag <- 0; LRV.ST <- 1
    }
    if (missing(LRV.kernel) || !is.character(LRV.kernel)) stop("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS' when structure is 'ts'")
    if (LRV.kernel=='uniform') weight.fn <- uniform.fn else if (LRV.kernel=='Bartlett') weight.fn <- Bartlett.fn else if (LRV.kernel=='QS') weight.fn <- QS.fn else stop(sprintf("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS'; not %s",LRV.kernel))
    # Compute gni() matrix
    gni.mat <- mZ*array(data=Itilde(-Lambda(tau=tau, theta=theta.hat, Y=Y, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi)/h)-tau,dim=dim(mZ))
    #
    if (is.na(LRV.ST)) { # Set ST automatically
      rho.hats <- sigma.hats <- rep(NA,dim(mZ)[2])
      for (a in 1:length(rho.hats)) {
        rho.hats[a] <- sum(gni.mat[1:(n-1),a]*gni.mat[2:n,a]) / sum(gni.mat[1:(n-1),a]^2)
        sigma.hats[a] <- suppressWarnings(sqrt(var(gni.mat[,a]) * (1-rho.hats[a]^2)))
      }
      if (any(c(is.nan(c(rho.hats,sigma.hats)),is.na(c(rho.hats,sigma.hats))))) {
        LRV.ST <- n^(1/5); if (LRV.kernel=='uniform' || LRV.kernel=='Bartlett') LRV.ST <- n^(1/3)
        warning(sprintf("AR method from Andrews (1991) for selecting S_T returned NA or NaN values; using S_T=%g",LRV.ST))
      } else if (LRV.kernel=='uniform') {
        LRV.ST <- uniform.ST.fn(alpha2=alpha2.fn(rhos=rho.hats,sigmas=sigma.hats),n=n)
        # stop("Need to provide a number for LRV.ST if LRV.kernel is 'uniform'")
      } else if (LRV.kernel=='Bartlett') {
        LRV.ST <- Bartlett.ST.fn(alpha1=alpha1.fn(rhos=rho.hats,sigmas=sigma.hats),n)
      } else if (LRV.kernel=='QS') {
        LRV.ST <- QS.ST.fn(alpha2=alpha2.fn(rhos=rho.hats,sigmas=sigma.hats),n)
      } else stop("Uncaught case.")
    }
    if (LRV.kernel!='QS') LRV.lag <- floor(LRV.ST)
    #
    tmpsum <- array(0,dim=rep(dim(mZ)[2],2))
    for (i in 1:n) {
      if (LRV.kernel=='QS') krange <- 1:n else krange <- max(1,i-LRV.lag):min(n,i+LRV.lag)
      for (k in krange) {
        tmpsum <- tmpsum + 
          weight.fn((i-k)/LRV.ST) * 
          (matrix(gni.mat[i,],ncol=1) %*% matrix(gni.mat[k,],nrow=1))
      }
    }
    return(tmpsum/(n-length(theta.hat))) #denominator adjustment per Andrews (1991) eqn (2.5)
  } else if (structure=='cluster') {
    stop("Not yet implemented: clustered covariance estimation")
  } else stop(sprintf("Argument structure must be either 'iid' or 'ts' or 'cluster' but its value is %s",structure))
}





























