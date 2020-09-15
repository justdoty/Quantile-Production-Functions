# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
set.seed(123456)
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/

#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source("gmmq.R")
source("ivqr_gmm.R")
#Required for 1st step estimation
require(quantreg)
#Required for QGMM estimation variations
require(GenSA)
require(pracma)
require(dplyr)
#Required for IVQR.BW (for two-step estimation)
require(MASS)
#Required for various codes to check ACF and LP estimates and bootstrapping
require(prodest)
QLP <- function(tau, va, state, free, proxy, idvar, timevar, h=0, b.init=NULL, R=20){
  #Make all data arguments into matrices
  va <- as.matrix(va)
  state <- as.matrix(state)
  free <- as.matrix(free)
  proxy <- as.matrix(proxy)
  idvar <- as.matrix(idvar)
  timevar<- as.matrix(timevar)
  snum <- ncol(state)
  fnum <- ncol(free)

  polyframe <- data.frame(state, proxy) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(free, mod, state^2, proxy^2)

  #Generate state lags
  lagstate <- state
  for (lag in 1:ncol(state)){
    lagstate[,lag] <- lagPanel(idvar=idvar, timevar=timevar, value=state[,lag])
  }
  #Generate free lags
  lagfree <- free
  for (lag in 1:ncol(free)){
    lagfree[,lag] <- lagPanel(idvar=idvar, timevar=timevar, value=free[,lag])
  }
  #Generate the matrix of data
  data <- suppressWarnings(as.matrix(data.frame(state=state, lagState=lagstate, free=free, lagFree=lagfree, va=va, idvar=idvar, timevar=timevar, regvars=regvars)))

  noboot <- finalQLP(tau=tau, h=h, ind=TRUE, data=data, fnum=fnum, snum=snum, b.init=b.init, boot=FALSE, gbar=TRUE)
  beta <- noboot[[1]]
  truegbar <- noboot[[2]]
  boot.indices <- block.boot.resample(idvar, R)
  boot.betas <- matrix(NA, R, (fnum+snum))

  for (i in 1:R){
    print(i)
    boot.betas[i,] <- finalQLP(tau=tau, h=h, ind=boot.indices[[i]], data=data, fnum=fnum, snum=snum, b.init=b.init, boot=TRUE, gbar=truegbar)
  }
  #Bootstrap estimates
  boot.beta <- apply(boot.betas, 2, mean)
  #Bootstrap standard errors
  boot.errors <- apply(boot.betas, 2, sd)
  #Returns to Scale:
  boot.scale <- mean(apply(boot.betas, 1, sum))
  #Returns to Scale standard error
  boot.scale.error <- sd(apply(boot.betas, 1, sum))
  #This returns the non-bootstrapped estimates, the matrix of bootstrapped estimates
  #The bootstrapped standard errors and the bootstrapped CI
  return(list(boot.betas, beta))
}
  #Function to estimate and to bootstrap QLP
finalQLP <- function(tau, h, ind, data, fnum, snum, b.init, boot, gbar){
  if (sum(as.numeric(ind))==length(ind)){ #if the ind variable is not always TRUE
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  #change the index according to bootstrapped indices
  data <- data[ind,] 
  first.stage <- rq(data[,'va', drop = FALSE] ~ data[, grepl('regvars', colnames(data)), drop = FALSE], na.action = na.exclude, tau=tau)
  free <- data[, grepl('free', colnames(data), ignore.case=FALSE), drop = FALSE]
  phi <- fitted(first.stage)
  beta.free <-  as.numeric(coef(first.stage)[2:(1+fnum)])
  #If not specified, starting points are the first stage+normal noise
  if (is.null(b.init)){
    b.init <- coef(first.stage)[(2 + fnum):(1 + fnum + snum)]
  } 
  #Clean Phi from the effects of free variables
  phi <- phi - (free%*%beta.free)
  newtime <- data[,'timevar', drop = FALSE]
  rownames(phi) <- NULL
  rownames(newtime) <- NULL
  #Lag Fitted Values  
  lag.phi <- lagPanel(value=phi, idvar=newid, timevar=newtime)
  #Clean the output from the effect of free variables
  va <- data[,'va', drop = FALSE] - (free%*%beta.free)
  state <- data[, grepl('state', colnames(data)), drop = FALSE]
  lagState <- data[, grepl('lagState', colnames(data)), drop = FALSE]
  lagFree <- data[, grepl('lagFree', colnames(data)), drop = FALSE]
  tmp.data <- na.omit(data.frame(state, lagState, lagFree, phi, lag.phi, va))
  lagFree <- tmp.data[, grepl('lagFree', colnames(tmp.data)), drop = FALSE]
  lagState <- tmp.data[, grepl('lagState', colnames(tmp.data)), drop = FALSE]
  #Instruments
  Z <- as.matrix(tmp.data$state)
  #Use W.init as an initial weighting matrix to get consistent estimates of theta
  # W.init <- solve(tau*(1-tau)*crossprod(Z)/nrow(Z))
  theta.hat <- as.numeric(ivqr.gmm(tau=tau, mX=tmp.data$state, mlX=tmp.data$lagState, mZ=Z, 
      fitphi=tmp.data$phi, fitlagphi=tmp.data$lag.phi, h=h, max.time=1, upper=1, lower=0, structure='iid', LRV.kernel='uniform', Lambda=Lambda, theta.init=b.init)$theta)
  #At the initial consistent estimate, compute the weighting matrix
  # W <- solve(LRV.est.fn(tau=tau, mX=tmp.data$state, mlX=tmp.data$lagState, 
  #   mZ=Z, fitphi=tmp.data$phi, fitlagphi=tmp.data$lag.phi, Lambda=Lambda, theta=theta.hat, 
  #   Itilde=Itilde.KS17, h=h, structure='iid', LRV.kernel='uniform'))
  #Solve using the above weighting matrix
  # try.state <- GenSA(par=theta.hat, fn=goQLP, mZ=Z, mW=W, mX=tmp.data$state, mlX=tmp.data$lagState, 
  # fitphi=tmp.data$phi, fitlagphi=tmp.data$lag.phi, tau=tau, h=h, gbar=gbar, lower=0, upper=1, control=list(max.time=5))
  # beta.state <- as.numeric(try.state$par)
  beta.state <- theta.hat

  if (gbar==TRUE){
    gbartrue <- g.bar(theta=beta.state, mX=tmp.data$state, mlX=tmp.data$lagState, mZ=Z, 
      fitphi=tmp.data$phi, fitlagphi=tmp.data$lag.phi, h=h, tau=tau)
    return(list(c(beta.state, beta.free), gbartrue))
  } else {
    return(c(beta.state, beta.free))
  }
}

#Smoothing Kernels
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
#Add bandwidth as second argument to G() and G'() functions.
  Gfn <- function(v,h){      
    Itilde.KS17(v/h)    
    }
  Gpfn <- function(v,h){      
    Itilde.deriv.KS17(v/h)    
    }

g.bar <- function(theta, mX, mlX, mZ, fitphi, fitlagphi, h, tau){
  theta <- as.matrix(as.numeric(theta))
  Omega <- fitphi-mX%*%theta
  Omega_lag <- fitlagphi-mlX%*%theta
  conc <- rq(Omega~Omega_lag, tau=0.5)
  rho <- as.numeric(coef(conc))[2]
  residconc <- resid(conc)
  beta0 <- quantile(residconc, tau)
  concparam <- as.numeric(c(beta0, rho))
  xifit <- Omega-cbind(1, Omega_lag)%*%concparam
  g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-xifit, h)-tau), 1, ncol(mZ))))
  return(g.bar)
}

#QLP Objective Function
goQLP <- function(tau, h, theta, mZ, mW, mX, mlX, fitphi, fitlagphi, gbar){
  if (gbar==TRUE){
    g.bar <- g.bar(theta, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau)
    crit <- nrow(mZ)*t(g.bar)%*%mW%*%g.bar
    return(crit)
  } else {
    g.bar.boot <- g.bar(theta, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau)-gbar
    return(nrow(mZ)*t(g.bar.boot)%*%mW%*%g.bar.boot)
  }
} 
Lambda <- function(theta, mX, mlX, fitphi, fitlagphi, tau){
  theta <- as.matrix(as.numeric(theta))
  Omega <- fitphi-mX%*%theta
  Omega_lag <- fitlagphi-mlX%*%theta
  conc <- rq(Omega~Omega_lag, tau=0.5)
  rho <- as.numeric(coef(conc))[2]
  residconc <- resid(conc)
  beta0 <- quantile(residconc, tau)
  concparam <- as.numeric(c(beta0, rho))
  xifit <- Omega-cbind(1, Omega_lag)%*%concparam
  return(xifit)
}  





