source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/ivqr_gmm.R')
require(dplyr)
###################################################################################################
#This file contains auxillary files for data subsetting and other miscellaneous tasks
#####################################################################################################

#####################################################################################################
#This function lags the data and returns both the lagged variables as well as its contemporaneous values
#####################################################################################################
lagdata <- function(idvar, X){
	condata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-1)
	lagdata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-n())
	data <- cbind(data.frame(condata), data.frame(lagdata[,-1]))
	return(data)
}
########################################################################################################
#boot resampling on IDs: bootstrapping on individuals (prodest.R)
#######################################################################################################
block.boot.resample <- function( idvar, R ){
  unique.ids <- unique(idvar) # find the unique values of panels in order to reshape the data
  panel.time.indices <- apply(unique.ids, 1, function(x) {return(list(which(idvar == x)))}) # find the time indices for each panel
  seq.indices <- 1:length(unique.ids) # the panel.time.indices list is indexed with sequential numbers: we mimic it
  boot.panel.id <- replicate(R, sample(seq.indices, replace = TRUE)) # generate a matrix of new IDs - R times
  new.indices <- list() # generate the matrix of the new indices
  ind <- 1:length(unique.ids)
  for (r in 1:R){ # for each boot rep we generate a vector of indices with rownames equal to a new - and fake - ID
    new.indices[[r]] <- cbind(unlist(mapply(function(x,y) {
      names(panel.time.indices[[x]][[1]]) <- rep(y,length(panel.time.indices[[x]][[1]]))
      return(list(panel.time.indices[[x]][[1]]))
    }, boot.panel.id[,r], ind))) # return a fake ID (sequential number) as row name and the index referring to the true ID
  }
  return(new.indices)
}

##################################################################################################
# Smoothing function Itilde() and its derviative Itilde'(); Vector input supported
##################################################################################################
Itilde.KS17 <- function(u) {  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0)) }
Itilde.deriv.KS17 <- function(u) { ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0) }
Gfn <- function(v,h){      
      Itilde.KS17(v/h)    
    }
Gpfn <- function(v,h){      
  Itilde.deriv.KS17(v/h)    
}

################################################################################################
#Kernel functions from Andrews (1991) eqn (2.7)
################################################################################################
QS.fn <- function(x) ifelse(x==0,1,(25/(12*pi^2*x^2))*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)))
Bartlett.fn <- function(x) ifelse(abs(x)>=1,0,1-abs(x))
uniform.fn <- function(x) ifelse(abs(x)>=1,0,1) #a.k.a. "Truncated" (Andrews 1991)
################################################################################################
# Optimal bandwidths from (5.2) (or 6.2) in Andrews (1991)
################################################################################################
QS.ST.fn <- function(alpha2,n) 1.3221*(alpha2*n)^(1/5)
Bartlett.ST.fn <- function(alpha1,n) 1.1447*(alpha1*n)^(1/3)
uniform.ST.fn <- function(alpha2,n) 0.6611*(alpha2*n)^(1/5) #footnote 5
################################################################################################
# (6.4) in Andrews (1991)
alpha2.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/(1-rhos)^8) / sum(ws*sigmas^4/(1-rhos)^4)
alpha1.fn <- function(rhos,sigmas,ws=1) sum(ws*4*rhos^2*sigmas^4/((1-rhos)^6*(1+rhos)^2)) / sum(ws*sigmas^4/(1-rhos)^4)

######################################################################################################
#Function that gives an estimate of the weighting matrix
######################################################################################################
LRV.est.fn <- function(tau, mX, mlX, mZ, fitphi, fitlagphi, theta, Lambda, Itilde, h, structure=c('iid','ts','cluster'), cluster.X.col, LRV.kernel=c('QS','Bartlett','uniform'), LRV.ST=NA, VERBOSE=FALSE) {
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
    gni.mat <- mZ*array(data=Itilde(-Lambda(theta=theta, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)/h)-tau,dim=dim(mZ))
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
    return(tmpsum/(n-length(theta))) #denominator adjustment per Andrews (1991) eqn (2.5)
  } else if (structure=='cluster') {
    stop("Not yet implemented: clustered covariance estimation")
  } else stop(sprintf("Argument structure must be either 'iid' or 'ts' or 'cluster' but its value is %s",structure))
}
###########################################################################################################
###########################################################################################################
###########################################################################################################