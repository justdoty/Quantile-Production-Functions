#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('PFQR/FUN/gmmq_aux.R')
#Required for QGMM estimation variations
require(pracma)
require(dplyr)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
LP <- function(idvar, timevar, Y, K, L, proxy, dZ, binit=NULL, R=20, tfptau){
  seed <- 123456
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
  trueboot <- finalQLP(ind=TRUE, data=data, binit=binit, gbar=0, seed=seed, tfptau=tfptau)
  #Elasticities from LP
  betahat <- trueboot$betahat
  #TFP quantiles
  QTFPhat <- trueboot$QTFPhat
  #Sample Moments used for recentering
  gbartrue <- trueboot$gbar
  #Initialize bootstrap#############################################################################
  #Indices used for resampling firm ID's
  bootind <- block.boot.resample(idvar, R, seed)
  #Bootstrapped elasticities for QLP
  betaboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrapped quantiles of TFP
  QTFPboot <- matrix(0, nrow=R, ncol=length(tfptau))
  #Bootstrap Procedure: finalQLP now computes the beta estimates where the sample moments
  #from GMM are recentered by truegbar, the sample moments evaluated at the true data
  for (i in 1:R){
    print(i)
    seed <- seed+i
    boot <- finalQLP(ind=bootind[[i]], data=data, binit=binit, gbar=gbartrue, seed=seed, tfptau=tfptau)
    betaboot[i,] <- boot$betahat
    QTFPboot[i,] <- boot$QTFPhat
  }
  return(list(betahat=betahat, QTFPhat=QTFPhat,  betaboot=betaboot, QTFPboot=QTFPboot))
}

###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(ind, data, binit, gbar, seed, tfptau){
  set.seed(seed)
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
  #Dataframe for regressors in first stage polynomial
  regvars  <- data.frame(reg1=data$K, reg2=data$L, reg3=data$proxy, reg4=data$proxy*data$K, reg5=data$K^2, reg6=data$proxy^2, reg7=data$proxy*(data$K^2), reg8=(data$proxy^2)*data$K, reg9=data$K^3, reg10=data$proxy^3)
  firststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  #Linear regression estimates (used for starting values)
  LM <- lm(data$Y~data$K+data$L)
  #Labor Estimate for LP
  LPLabor <-  as.numeric(coef(firststage)[3])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LPLabor
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  #Output net of labor
  mY <- as.matrix(newdata$Ycon-newdata$Lcon*LPLabor)
  #Contemporary State Variables
  mX <- as.matrix(newdata$Kcon)
  #Lagged State Variables
  mlX <- as.matrix(newdata$Klag)
  #Contemporary phi estimates
  fitphi <- as.matrix(newdata$phicon)
  #Lagged phi estimates
  fitlagphi <- as.matrix(newdata$philag)
  #Instruments
  mZ <- cbind(as.matrix(newdata$Kcon), as.matrix(newdata$Klag), as.matrix(newdata$Llag), as.matrix(newdata$Pxlag))
  #If not specified, starting point is the first stage estimates
  if (is.null(binit)){
    binit <- as.numeric(coef(LM)[2])

  }
  if (ncol(mZ)>ncol(mX)){
    soln <- optim(par=binit, fn=function(b) LPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    xhat <- soln$par
    betahat <- c(xhat[1], LPLabor)
    #Value of Moments
    gbar <- gbar(b=soln$par, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar)
    #TFP estimates (in logs)
    TFP <- data$Y-cbind(data$K, data$L)%*%betahat
    #Quantiles of log TFP
    QTFPhat <- as.numeric(quantile(TFP, tfptau))
    return(list(betahat=betahat, QTFPhat=QTFPhat, gbar=gbar))
    
  } else if (ncol(mZ)==ncol(mX)){
    soln <- optim(par=binit, fn=function(b) LPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    xhat <- soln$par
    betahat <- c(xhat[1], LPLabor)
    #Value of Moments (Set at 0 since no recentering)
    gbar <- 0
    #TFP estimates (in logs)
    TFP <- data$Y-cbind(data$K, data$L)%*%betahat
    #Quantiles of log TFP
    QTFPhat <- as.numeric(quantile(TFP, tfptau))
    return(list(betahat=betahat, QTFPhat=QTFPhat, gbar=gbar))
  }
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- mY-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B+I(B^2)+I(B^3))
  step1param <- as.numeric(coef(step1))
  xifit <- A-cbind(1, B, B^2, B^3)%*%step1param
  return(xifit)
} 
##################################################################################
#This function calculates the sample moment for a given b
##################################################################################
gbar <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi, gbar){
  xifit <- Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
  gbar <- as.matrix(colMeans(gni))
  return(gbar)
}
################################################################################
#LP Objective Function
################################################################################
LPobj <- function(b, mY, mX, mZ, mlX, fitphi, fitlagphi, gbar){
  xifit <- Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
  gnic <- colMeans(gni)
  go <- mean(gnic^2)
  return(go)
} 
##################################################################################
##################################################################################
##################################################################################














