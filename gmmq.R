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
gmmq <- function(tau, X, Z, dB, 
                 Lambda=function(x,b) x[,1]-x[,-1]%*%b, 
                 Lambda.derivative=function(y,x,b) -x[,-1], 
                 h=0, VERBOSE=FALSE, RETURN.Z=FALSE, b.init=0) {
  if (!require(pracma)) stop('Please install (from CRAN) the R package pracma to run this code.')
  # Check arguments: errors, warnings, defaults
  if (missing(tau) || missing(X)) stop('Must supply tau (scalar), dB (scalar), X (matrix)')
  if (!is.numeric(tau)|| !is.numeric(X) || (!is.null(Z) && !is.numeric(Z)) || !is.numeric(h)) {
    stop('All arguments besides Lambda must be numeric.')
  }
  if (!is.null(Z) && !is.array(Z)) Z <- matrix(Z, ncol=1)
  if (missing(dB)) stop('Enter the number of parameters to be estimated')
  if (!is.numeric(dB)) stop("dB must be numeric")
  if (b.init!=0 && length(b.init)!=dB) stop(sprintf("Length of b.init must equal dB, which is %d", dB))
  if (!is.function(Lambda)) stop('Lambda must be a function.')
  if (!is.null(Lambda.derivative) && !is.function(Lambda.derivative)) stop('Lambda.derivative must be a function (or NULL).')
  if (!is.array(X) || (!is.null(Z) && !is.array(Z))) stop('X, and Z (if not NULL) must be arrays/matrices.')
  if (length(tau)!=1 || length(dB)!=1 || length(h)!=1) stop('tau, dB, and h must all be scalars: length(dB)==1, etc.')
  if (dim(X)[1]!=dim(Z)[1]) stop('X and Z must have same number of rows.')
  if (tau<=0 || tau>=1) stop('Quantile index tau (first argument) must be strictly between 0 and 1.')
  if (h<0) stop('Initial bandwidth h must be nonnegative.')
  # Check for NA, remove rows
  NAany <- which(apply(X=Z,MARGIN=1,FUN=function(row)any(is.na(row))))
  if (length(NAany)>0) {
    cat(sprintf("Removing %d observations due to NA, X, or Z\n",length(NAany)))
    Z <- matrix(Z[-NAany,], ncol=dim(Z)[2])
    X <- matrix(X[-NAany,], ncol=dim(X)[2])
  }
  n <- dim(X)[1]
  #
  if (dB>dim(Z)[2]) stop('Parameters not identified: must have dB<=dim(X)[2]+dim(Z)[2]') 
 
  # Set smoothing function Itilde() and its derviative Itilde'()
  # (Vector input supported)
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17

  # For solver (to find solution to estimating equations)
  # YX <- cbind(Y,X)
  obj.fn <- function(b,h) { #objective function
    # L <- matrix(apply(cbind(Y,X),1,function(YXobs)Lambda(YXobs[1:dY],YXobs[-(1:dY)],b)), ncol=1)
    # L <- matrix(Lambda(y=Y,x=X,b=b), ncol=1)
    L <- Lambda(x=X,b=b)
    return(colMeans(Z*array(data=Itilde(-L/h)-tau,dim=dim(Z))))
  }
  if (!is.null(Lambda.derivative)) {
    jac.fn <- function(b,h) { #Jacobian function
      # L <- matrix(apply(cbind(Y,X),1,function(YXobs)Lambda(YXobs[1:dY],YXobs[-(1:dY)],b)), ncol=1)
      # L <- matrix(Lambda(y=Y,x=X,b=b), ncol=1)
      L <- Lambda(x=X,b=b)
      # Lp <- t(apply(YX,1,function(YXobs)Lambda.derivative(YXobs[1:dY],YXobs[-(1:dY)],b)))
      Lp <- Lambda.derivative(x=X,b=b)
      # if (VERBOSE) cat(sprintf("dim(L)=%s;dim(Lp)=%s;dB=%d\n",paste0(sprintf("%d",dim(L)),collapse="-by-"),paste0(sprintf("%d",dim(Lp)),collapse="-by-"),dB))
      return((t(Z) %*% (Lp*array(data=Itilde.deriv(-L/h),dim=dim(Z))*(-1/h)))/n)
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
  b.best <- matrix(data=b.init,nrow=dB,ncol=1)
  MAXITER <- 400;  H.MAX <- 1e10 #.Machine$double.xmax/1e3
  # Loop: h<=h.lo<h.cur<h.hi, where h.hi can be solved by newtonsys but h.lo cannot
  if (VERBOSE) cat(sprintf('n=%d, tau=%g, dB=%d, desired h=%g\n',n,tau,dB,h))
  last.error <- NULL
  while (h.cur<=H.MAX && h.cur>.Machine$double.eps*1e2 && 
         (is.na(h.hi) || (h.hi/h.lo)>1.4)) {
    if (VERBOSE) cat(sprintf("h.cur=%g: ",h.cur))
    if (is.null(Lambda.derivative)) J <- NULL else J <- function(b)jac.fn(b,h.cur)
  	soln <- tryCatch(newtonsys(Ffun=function(b)obj.fn(b,h.cur),Jfun=J,
  	                            x0=b.best,maxiter=MAXITER),
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
#   	  } else if (h.lo==0 || h.lo==h) { #maybe with better x0 in newtonsys, a solution can be found
#   		  h.cur <- h.init
#   		} else { #try geometric mean of h.lo and h.hi
#   		  h.cur <- sqrt(h.lo*h.hi)
#   		}
  	} 
  }
  if (h.cur>H.MAX) {
    warning("(above) error from pracma::newtonsys")
    # stop(last.error)
  }
  if (RETURN.Z) return(list(b=b.best,h=h.hi,Z=Z)) else return(list(b=b.best,h=h.hi))
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

LRV.est.fn <- function(tau,X,Z,Lambda,beta.hat,Itilde,h,structure=c('iid','ts','cluster'),cluster.X.col,LRV.kernel=c('QS','Bartlett','uniform'),LRV.ST=NA,VERBOSE=FALSE) {
  # if (missing(structure) || !is.character(structure)) stop("Argument structure must be 'iid' or 'ts' or 'cluster'")
  structure <- match.arg(structure)
  LRV.kernel <- match.arg(LRV.kernel)
  n <- dim(Z)[1]
  if (structure %in% c('iid','ts')) {
    if (structure=='iid') {
      LRV.kernel <- 'uniform'; LRV.lag <- 0; LRV.ST <- 1
    }
    if (missing(LRV.kernel) || !is.character(LRV.kernel)) stop("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS' when structure is 'ts'")
    if (LRV.kernel=='uniform') weight.fn <- uniform.fn else if (LRV.kernel=='Bartlett') weight.fn <- Bartlett.fn else if (LRV.kernel=='QS') weight.fn <- QS.fn else stop(sprintf("LRV.kernel must be 'uniform' or 'Bartlett' or 'QS'; not %s",LRV.kernel))
    # Compute gni() matrix
    gni.mat <- Z*array(data=Itilde(-Lambda(x=X,b=beta.hat)/h)-tau,dim=dim(Z))
    #
    if (is.na(LRV.ST)) { # Set ST automatically
      rho.hats <- sigma.hats <- rep(NA,dim(Z)[2])
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
    tmpsum <- array(0,dim=rep(dim(Z)[2],2))
    for (i in 1:n) {
      if (LRV.kernel=='QS') krange <- 1:n else krange <- max(1,i-LRV.lag):min(n,i+LRV.lag)
      for (k in krange) {
        tmpsum <- tmpsum + 
          weight.fn((i-k)/LRV.ST) * 
          (matrix(gni.mat[i,],ncol=1) %*% matrix(gni.mat[k,],nrow=1))
      }
    }
    return(tmpsum/(n-length(beta.hat))) #denominator adjustment per Andrews (1991) eqn (2.5)
  } else if (structure=='cluster') {
    stop("Not yet implemented: clustered covariance estimation")
  } else stop(sprintf("Argument structure must be either 'iid' or 'ts' or 'cluster' but its value is %s",structure))
}






























