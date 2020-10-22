# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('gmmq_aux.R')
#Required for 1st step estimation
require(quantreg)
#Required for QGMM estimation variations
require(GenSA)
require(pracma)
require(dplyr)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
QLP <- function(tau, idvar, timevar, Y, K, L, proxy, h=0, b.init=NULL, R=20){
  seed <- 123456
  set.seed(seed)
  #Make all data arguments into matrices
  idvar <- as.matrix(idvar)
  timevar <- as.matrix(timevar)
  Y <- as.matrix(Y)
  K <- as.matrix(K)
  L <- as.matrix(L)
  proxy <- as.matrix(proxy)
  #####################
  data <- data.frame(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, proxy=proxy)
  #This function computes the "true" beta and sample moments evaluated at the "true" parameters
  #using the "true" data used for recentering the moments in the bootstrap
  #Here ind denotes the index that does tells finalQLP not resample the original data
  trueboot <- finalQLP(tau=tau, h=h, ind=TRUE, data=data, b.init=b.init, gbar=0)
  #True parameters
  beta <- trueboot$beta
  #True sample moments
  gbar <- trueboot$gbartrue
  #Initialize bootstrap
  boot.indices <- block.boot.resample(idvar, R)
  boot.betas <- matrix(0, nrow=R, ncol=2)
  #Bootstrap Procedure: finalQLP now computes the beta estimates where the sample moments
  #from GMM are recentered by truegbar, the sample moments evaluated at the true data
  for (i in 1:R){
    print(i)
    set.seed(seed+i)
    boot.betas[i,] <- finalQLP(tau=tau, h=h, ind=boot.indices[[i]], data=data, b.init=b.init, gbar=gbar)$beta
  }
  #Only return the boostrapped estimates and the estimates at the true data
  return(list(boot.betas, beta))
}

###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(tau, h, ind, data, b.init, gbar){
  ##########################################################################
  #Bootstrap Component
  ###########################################################################
  if (sum(as.numeric(ind))==length(ind)){ 
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  data <- data[ind,] 
  #######################################################################
  #First Stage of LP
  #######################################################################
  #Create a polynomail data frame for 1st stage estimation
  regvars <- data.frame(reg1=data$L, reg2=data$K, reg3=data$K*data$proxy, reg4=data$K^2, reg5=data$proxy^2, reg6=data$K^2*data$proxy^2)
  firststage <- rq(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]), tau=tau)
  phi0 <- as.numeric(coef(first.stage)[1])
  LP_Labor <-  as.numeric(coef(first.stage)[2])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor-phi0
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  #Output net of labor
  mY <- newdata$Ycon-newdata$Lcon*LP_Labor
  #Instruments
  Z <- cbind(1, as.matrix(newdata$Kcon))
  #Use W.init as an initial weighting matrix to get consistent estimates of theta
  Winit <- solve(tau*(1-tau)*crossprod(Z)/nrow(Z))
  #If not specified, starting point is the first stage estimates
  if (is.null(b.init)){
    b.init <- as.numeric(coef(first.stage)[3])
  } 
  #Solve using the above weighting matrix
  if (ncol(Z)>ncol(newdata$Kcon)){
    kinit <- GenSA(par=b.init, fn=QLPgo, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=Winit, gbartrue=gbar, h=h, tau=tau, lower=0, upper=1, control=list(max.time=5))$par
    W <- solve(LRV.est.fn(tau=tau, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, 
    mZ=Z, fitphi=newdata$phicon, fitlagphi=newdata$philag,, Lambda=Lambda, theta=kinit, 
    Itilde=Itilde.KS17, h=h, structure='iid', LRV.kernel='uniform'))
    khat <- GenSA(par=kinit, fn=QLPgo, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=W, gbartrue=gbar, h=h, tau=tau, lower=0, upper=1, control=list(max.time=5))$par
    gbartrue <- gbar(theta=khat, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fithpi=newdata$phicon,
      lagphi=newdata$philag, h=h, tau=tau)
    return(list(beta=c(khat, LP_Labor), gbartrue=gbartrue))
    
  } else if (ncol(Z)==ncol(newdata$Kcon)){
    kinit <- GenSA(par=b.init, fn=QLPgo, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=Winit, h=h, gbartrue=gbar, tau=tau, lower=0, upper=1, control=list(max.time=5))$par
    W <- solve(LRV.est.fn(tau=tau, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, 
    mZ=Z, fitphi=newdata$phicon, fitlagphi=newdata$philag,, Lambda=Lambda, theta=kinit, 
    Itilde=Itilde.KS17, h=h, structure='iid', LRV.kernel='uniform'))
    khat <- GenSA(par=kinit, fn=QLPgo, mY=mY, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=W, h=h, gbartrue=gbar, tau=tau, lower=0, upper=1, control=list(max.time=5))$par
    gbartrue <- 0
    return(list(beta=c(khat, LP_Labor), gbar=gbartrue))
  }
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(theta, mY, mX, mlX, fitphi, fitlagphi, tau){
  theta <- as.matrix(as.numeric(theta))
  A <- mY-mX%*%theta[1:ncol(mX)]
  B <- fitlagphi-mX%*%theta[1:ncol(mX)]
  step1 <- rq(A~B+B^2+B^3, tau=0.5)
  step1param <- as.numeric(coef(step1))
  residconc <- A-cbind(B, B^2, B^3)%*%step1param[-1]
  Finv <- quantile(residconc, tau)
  xifit <- residconc-Finv
  return(xifit)
} 
##################################################################################
#This function calculates the smoothed sample moment for a given theta
##################################################################################
gbar <- function(theta, mY, mX, mlX, mZ, fitphi, fitlagphi, h, tau){
  xifit <- Lambda(theta=theta, mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, h=h, tau=tau)
  gbar <- as.matrix(colMeans(mZ*repmat((Gfn(-xifit, h)-tau), 1, ncol(mZ))))
  return(gbar)
}
################################################################################
#QLP Objective Function
################################################################################
QLPgo <- function(theta, mY, mX, mlX, fitphi, fitlagphi, mW, gbartrue, h, tau){
  gbar <- gbar(theta, mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau)-gbartrue
  go <- nrow(mZ)*t(gbar)%*%mW%*%gbar
  return(go)
} 
##################################################################################
##################################################################################
##################################################################################




