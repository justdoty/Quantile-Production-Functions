#We use the bootstrap code in prodest.R, an R package for production function estimation written by Gabrielle Rovigatti
source('PFQR/FUN/Aux_Fun.R')
require(quantreg)
require(dplyr)
require(pracma)
require(GenSA)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
QLP_Boot <- function(tau, idvar, timevar, Y, K, L, proxy, dZ, binit=NULL, R=20){
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
  trueboot <- finalQLP(tau=tau, ind=TRUE, data=data, binit=binit, seed=seed)
  #Estimates from the un-resampled###############################################################
  #Elasticities from QLP
  betahat <- trueboot$betahat
  LPhat <- trueboot$LPhat
  #Estimates from QR
  qrhat <- trueboot$qrhat
  #Difference between QLP and QR
  qdifhat <- trueboot$qdifhat
  #QTFP estimates
  QTFPhat <- trueboot$QTFPhat
  #TFP estimates
  TFPhat <- trueboot$TFPhat
  #Omega estimates
  omegahat <- trueboot$omegahat
  expost <- trueboot$expost
  #Initialize bootstrap#############################################################################
  #Indices used for resampling firm ID's
  bootind <- block.boot.resample(idvar, R, seed)
  #Bootstrapped elasticities for QLP
  betaboot <- matrix(0, nrow=R, ncol=dZ)
  LPboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrapped QR estimates
  qrboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrapped differences between QLP and QR
  qdifboot <- matrix(0, nrow=R, ncol=dZ)
  for (i in 1:R){
    print(i)
    seed <- seed+i
    boot <- finalQLP(tau=tau, ind=bootind[[i]], data=data, binit=binit, seed=seed)
    betaboot[i,] <- boot$betahat
    LPboot[i,] <- boot$LPhat
    qrboot[i,] <- boot$qrhat
    qdifboot[i,] <- boot$qdifhat
  }
  return(list(betahat=betahat, LPhat=LPhat, qrhat=qrhat, qdifhat=qdifhat, QTFPhat=QTFPhat, TFPhat=TFPhat, omegahat=omegahat,  betaboot=betaboot, LPboot=LPboot, 
    qrboot=qrboot, qdifboot=qdifboot, expost=expost))
}
###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(tau, ind, data, binit, seed){
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
  regvars <- data.frame(reg1=data$K, reg2=data$L, reg3=data$proxy, reg4=data$proxy*data$K, reg5=data$K^2, reg6=data$proxy^2, reg7=data$proxy*(data$K^2), reg8=(data$proxy^2)*data$K, reg9=data$K^3, reg10=data$proxy^3)
  #First stage estimates for LP
  LPfirststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  #Labor Estimate for LP
  LPLabor <-  as.numeric(coef(LPfirststage)[3])
  #QR estimates
  QR <- rq(data$Y~data$K+data$L, tau=tau)
  qrhat <- as.numeric(coef(QR)[-1])
  #Linear regression estimates (used for starting values)
  LM <- lm(data$Y~data$K+data$L)
  #Clean Phi from the effects of free variables
  LPphi <- fitted(LPfirststage)-as.matrix(data$L)%*%LPLabor
  #Ex-post shocks
  expost <- data$Y-fitted(LPfirststage)
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, LPphi))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "LPphicon", "Ylag", "Klag", "Llag", "Pxlag", "LPphilag")
  #LP Output net of labor
  LPmY <-  as.matrix(newdata$Ycon-newdata$Lcon*LPLabor)
  #LP Contemporary State Variables
  LPmX <- as.matrix(newdata$Kcon)
  #LP Lagged State Variables
  LPmlX <- as.matrix(newdata$Klag)
  #LP Contemporary phi estimates
  LPfitphi <- as.matrix(newdata$LPphicon)
  #LP Lagged phi estimates
  LPfitlagphi <- as.matrix(newdata$LPphilag)
  #LP Instruments (Exact Identification)
  LPmZ <- as.matrix(newdata$Kcon)
  #Starting values for LP estimates from OLS
  LPkinit <- as.numeric(coef(LM)[2])
  #LP estimates for Capital
  LPkhat <- optim(par=LPkinit, fn=function(b) LPobj(b, mY=LPmY, mX=LPmX, mlX=LPmlX, mZ=LPmZ, fitphi=LPfitphi, fitlagphi=LPfitlagphi), gr=NULL, method="L-BFGS-B", lower=0, upper=1)$par
  #Estimates of productivity from LP
  omegahat <- LPphi-data$K*LPkhat
  #Output net of productivity
  mY <- as.matrix(data$Y-omegahat)
  #State Variables
  mX <- cbind(data$K, data$L)
  #QLP Estimates for Capital from (QR) good starting values as well
  mom <- rq(mY~mX-1, tau=tau)
  betahat <- as.numeric(coef(mom))
  #LP Estimates
  LPhat <- c(LPkhat, LPLabor)
  #Difference between QLP and QR estimates
  qdifhat <- betahat-qrhat
  #QLP TFP estimates (in logs)
  QTFPhat <- data$Y-cbind(data$K, data$L)%*%betahat
  #LP TFP estimates (in logs)
  TFPhat <- data$Y-cbind(data$K, data$L)%*%LPhat
  return(list(betahat=betahat, LPhat=LPhat, qrhat=qrhat, qdifhat=qdifhat, QTFPhat=QTFPhat, TFPhat=TFPhat, omegahat=omegahat, expost=expost))
}
############################################################################################
#Functions for Estimating LP Coefficients
###########################################################################################
#Function that defines the residuals
LP_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B+I(B^2)+I(B^3))
  step1param <- as.numeric(coef(step1))
  wfit <- cbind(1, B, B^2, B^3)%*%step1param
  resid <- mY-mX%*%b[1:ncol(mX)]-wfit
  return(resid)
} 
#LP GMM objective function
LPobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi){
  resid <- LP_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  #Since capital is exogeneous in this model, we do not consider using instruments in estimation
  #Instead, simply use sum of squared errors to get a consistent estimate of capital elasticity
  xi <- crossprod(resid)
  # gni <- mZ*repmat(resid,1, ncol(mZ))
  # gnic <- colSums(gni^2)
  # xi <- sum(gnic)
  return(xi)
}



