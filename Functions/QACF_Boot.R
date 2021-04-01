#We use the bootstrap code in prodest.R, an R package for production function estimation written by Gabrielle Rovigatti
source('PFQR/FUN/Aux_Fun.R')
require(quantreg)
require(dplyr)
require(pracma)
require(GenSA)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQACF
###################################################################################
###################################################################################
QACF_Boot <- function(tau, idvar, timevar, Y, K, L, proxy, dZ, binit=NULL, R=20){
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
  #Here ind denotes the index that does tells finalACF not resample the original data
  trueboot <- finalQACF(tau=tau, ind=TRUE, data=data, binit=binit, seed=seed)
  #Estimates from the un-resampled###############################################################
  #Elasticities from QACF
  betahat <- trueboot$betahat
  ACFhat <- trueboot$ACFhat
  #Estimates from QR
  qrhat <- trueboot$qrhat
  #Difference between QACF and QR
  qdifhat <- trueboot$qdifhat
  #QTFP estimates
  QTFPhat <- trueboot$QTFPhat
  #TFP estimates
  TFPhat <- trueboot$TFPhat
  #Omega estimates
  omegahat <- trueboot$omegahat
  #Estimates of Ex-post shock
  expost <- trueboot$expost
  #Initialize bootstrap#############################################################################
  #Indices used for resampling firm ID's
  bootind <- block.boot.resample(idvar, R, seed)
  #Bootstrapped elasticities for QACF
  betaboot <- matrix(0, nrow=R, ncol=dZ)
  ACFboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrapped QR estimates
  qrboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrapped differences between QACF and QR
  qdifboot <- matrix(0, nrow=R, ncol=dZ)
  #Bootstrap Procedure:
  for (i in 1:R){
    print(i)
    seed <- seed+i
    boot <- finalQACF(tau=tau, ind=bootind[[i]], data=data, binit=binit, seed=seed)
    betaboot[i,] <- boot$betahat
    ACFboot[i,] <- boot$ACFhat
    qrboot[i,] <- boot$qrhat
    qdifboot[i,] <- boot$qdifhat
  }
  return(list(betahat=betahat, ACFhat=ACFhat, qrhat=qrhat, qdifhat=qdifhat, QTFPhat=QTFPhat, TFPhat=TFPhat, omegahat=omegahat,  betaboot=betaboot, ACFboot=ACFboot, 
    qrboot=qrboot, qdifboot=qdifboot, expost=expost))
}
###########################################################################
###########################################################################
#Function to estimate and to bootstrap QACF
###########################################################################
###########################################################################
finalQACF <- function(tau, ind, data, binit, seed){
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
  #First Stage of ACF
  #######################################################################
  #Dataframe for regressors in first stage polynomial
  regvars <- data.frame(reg1=data$K, reg2=data$L, reg3=data$proxy, reg4=data$K*data$L, reg5=data$K*data$proxy, 
    reg6=data$L*data$proxy, reg7=data$K^2, reg8=data$L^2, reg9=data$proxy^2)
  #First stage estimates for ACF
  ACFfirststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  #QR estimates
  QR <- rq(data$Y~data$K+data$L, tau=tau)
  qrhat <- as.numeric(coef(QR)[-1])
  #Linear regression estimates (used for starting values)
  LM <- lm(data$Y~data$K+data$L)
  #First stage fitted values
  phi0 <- as.numeric(coef(ACFfirststage)[1])
  phihat <- fitted(ACFfirststage)
  expost <- data$Y-phihat
  ACFphi <- phihat
  #Calculate Contempory and Lag Values for 2nd stage estimation
  lagdata1 <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, ACFphi))
  names(lagdata1) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "ACFphicon", "Ylag1", "Klag1", "Llag1", "Pxlag1", "ACFphilag1")
  lagdata2 <- lagdata(idvar=lagdata1$idvar, X=cbind(lagdata1$Ycon, lagdata1$Kcon, lagdata1$Lcon, 
    lagdata1$Pxcon, lagdata1$ACFphicon, lagdata1$Ylag, lagdata1$Klag, lagdata1$Llag, 
    lagdata1$Pxlag, lagdata1$ACFphilag))
  #Naming convention below is bad practice, but used to prevent the copy of lag1 variables in the new dataset
  names(lagdata2) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "ACFphicon", "Ylag1", "Klag1", "Llag1", "Pxlag1", "ACFphilag1",
    "NA", "NA", "NA", "NA", "NA", "Ylag2", "Klag2", "Llag2", "Pxlag2", "ACFphilag2")
  #ACF Output
  # ACFmY <-  as.matrix(lagdata1$Ycon)
  ACFmY <-  as.matrix(lagdata2$Ycon)
  #ACF Contemporary State Variables
  # ACFmX <- cbind(lagdata1$Kcon, lagdata1$Lcon)
  ACFmX <- cbind(lagdata2$Kcon, lagdata2$Lcon)
  #ACF Lagged State Variables
  # ACFmlX <- cbind(lagdata1$Klag1, lagdata1$Llag1)
  ACFmlX <- cbind(lagdata2$Klag1, lagdata2$Llag1)
  #ACF Contemporary phi estimates
  # ACFfitphi <- as.matrix(lagdata1$ACFphicon)
  ACFfitphi <- as.matrix(lagdata2$ACFphicon)
  #ACF Lagged phi estimates
  # ACFfitlagphi <- as.matrix(lagdata1$ACFphilag1)
  ACFfitlagphi <- as.matrix(lagdata2$ACFphilag1)
  #ACF Instruments 
  # ACFmZ <- cbind(lagdata1$Kcon, lagdata1$Llag1)
  ACFmZ <- cbind(1, lagdata2$Kcon, lagdata2$Llag1, lagdata2$Klag, lagdata2$Llag2)
  #Starting values for ACF estimates from OLS
  ACFinit <- as.numeric(coef(LM)[-1])
  #ACF estimates for Capital
  ACFhat <- optim(par=ACFinit, fn=function(b) ACFobj(b, mY=ACFmY, mX=ACFmX, mlX=ACFmlX, mZ=ACFmZ, fitphi=ACFfitphi, fitlagphi=ACFfitlagphi), gr=NULL, method="L-BFGS-B", lower=c(0,0), upper=c(1.5,1.5))$par
  #Estimates of productivity from ACF
  wfit <- phihat-cbind(data$K, data$L)%*%as.matrix(as.numeric(ACFhat))
  #Output net of productivity
  mY <- as.matrix(data$Y-wfit)
  #State Variables
  mX <- cbind(data$K, data$L)
  #QACF Estimates for Capital and Labor
  mom <- rq(mY~mX-1, tau=tau)
  betahat <- as.numeric(coef(mom))
  #Difference between QACF and QR estimates
  qdifhat <- betahat-qrhat
  #QACF TFP estimates (in logs)
  QTFPhat <- data$Y-cbind(data$K, data$L)%*%betahat
  #ACF TFP estimates (in logs)
  TFPhat <- data$Y-cbind(data$K, data$L)%*%as.matrix(as.numeric(ACFhat))
  #QACF Productivity Estimates (in logs and unweighted)
  omegahat <- wfit
  return(list(betahat=betahat, ACFhat=ACFhat, qrhat=qrhat, qdifhat=qdifhat, QTFPhat=QTFPhat, TFPhat=TFPhat, omegahat=omegahat, expost=expost))
} 
############################################################################################
#Functions for Estimating ACF Coefficients
###########################################################################################
#Function that defines the residuals
ACF_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  # step1 <- lm(A~B)
  step1 <- lm(A~B-1)
  step1param <- as.numeric(coef(step1))
  # xifit <- A-cbind(1, B)%*%step1param
  xifit <- A-B*step1param
  return(xifit)
} 
#ACF GMM objective function
ACFobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi){
  xifit <- ACF_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  mW <- solve(crossprod(mZ))/nrow(mZ)
  go <- t(crossprod(mZ, xifit))%*%mW%*%(crossprod(mZ, xifit))
  return(go)
}



