source('PFQR/FUN/Aux_Fun.R')
require(quantreg)
require(dplyr)
require(pracma)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQACF
###################################################################################
###################################################################################
QACF_Boot <- function(tau, idvar, timevar, Y, K, L, proxy, dZ, binit=NULL, R=20, XC, XB){
  seed <- 123456
  #Make all data arguments into matrices
  idvar <- as.matrix(idvar)
  timevar <- as.matrix(timevar)
  Y <- as.matrix(Y)
  K <- as.matrix(K)
  L <- as.matrix(L)
  proxy <- as.matrix(proxy)
  #####################
  data <- data.frame(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, proxy=proxy, XC=XC, XB=XB)
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
  #Productivity Differentials
  ACFPB <- trueboot$ACFPB
  DSPB <- trueboot$DSPB
  ACFPC <- trueboot$ACFPC
  DSPC <- trueboot$DSPC
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
  #Productivity Differentials
  ACFPBboot <- matrix(0, nrow=R, ncol=ncol(XB))
  DSPBboot <- matrix(0, nrow=R, ncol=ncol(XB))
  ACFPCboot <- matrix(0, nrow=R, ncol=ncol(XC))
  DSPCboot <- matrix(0, nrow=R, ncol=ncol(XC))
  #Bootstrap Procedure:
  for (i in 1:R){
    print(i)
    seed <- seed+i
    boot <- finalQACF(tau=tau, ind=bootind[[i]], data=data, binit=binit, seed=seed)
    betaboot[i,] <- boot$betahat
    ACFboot[i,] <- boot$ACFhat
    qrboot[i,] <- boot$qrhat
    qdifboot[i,] <- boot$qdifhat
    ACFPBboot[i,] <- boot$ACFPB
    DSPBboot[i,] <- boot$DSPB
    ACFPCboot[i,] <- boot$ACFPC
    DSPCboot[i,] <- boot$DSPC
  }
  return(list(betahat=betahat, ACFhat=ACFhat, qrhat=qrhat, qdifhat=qdifhat,  betaboot=betaboot, ACFboot=ACFboot, 
    qrboot=qrboot, qdifboot=qdifboot, ACFPB=ACFPB, DSPB=DSPB, ACFPC=ACFPC, DSPC=DSPC, ACFPBboot=ACFPBboot,
    DSPBboot=DSPBboot, ACFPCboot=ACFPCboot, DSPCboot=DSPCboot))
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
  ACFphi <- phihat
  #Calculate Contempory and Lag Values for 2nd stage estimation
  idcon <- duplicated(data$idvar)
  idcon1 <- duplicated(duplicated(data$idvar))
  idlag <- duplicated(data$idvar, fromLast=TRUE)
  idlag1 <- duplicated(duplicated(data$idvar, fromLast=TRUE))
  idlag2 <- duplicated(duplicated(data$idvar, fromLast=TRUE))
  #ACF Output
  ACFmY <-  data$Y[idcon]
  #ACF Contemporary State Variables
  ACFmX <- cbind(data$K[idcon], data$L[idcon])
  # ACFmX <- cbind(lagdata2$Kcon, lagdata2$Lcon)
  #ACF Lagged State Variables
  ACFmlX <- cbind(data$K[idlag], data$L[idlag])
  # ACFmlX <- cbind(lagdata2$Klag1, lagdata2$Llag1)
  #ACF Contemporary phi estimates
  ACFfitphi <- ACFphi[idcon]
  # ACFfitphi <- as.matrix(lagdata2$ACFphicon)
  #ACF Lagged phi estimates
  ACFfitlagphi <- ACFphi[idlag]
  # ACFfitlagphi <- as.matrix(lagdata2$ACFphilag1)
  #ACF Instruments 
  ACFmZ <- cbind(1, data$K[idcon], data$L[idlag])
  # ACFmZ <- cbind(1, lagdata2$Kcon, lagdata2$Llag1, lagdata2$Klag, lagdata2$Llag2)
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
  #Estimates of Total Factor Productivity from ACF
  ACFTFPhat <- data$Y-cbind(data$K, data$L)%*%as.matrix(ACFhat)
  #Estimates of Total Factor Productivity from DS
  DSTFPhat <- data$Y-cbind(data$K, data$L)%*%as.matrix(betahat)
  XC <- as.matrix(data[,grepl("XC", names(data))])
  XB <- as.matrix(data[,grepl("XB", names(data))])
  pdata <- data.frame(ACF=ACFTFPhat, DS=DSTFPhat, XC=XC, XB=XB) %>% na.omit()
  XC <- as.matrix(pdata[,grepl("XC", names(pdata))])
  XB <- as.matrix(pdata[,grepl("XB", names(pdata))])
  XCind <- apply(XC==0, 1, sum)<1
  XC <- XC[XCind,]
  #Estimates of Productivity Differentials from ACF
  ACFPB <- apply(XB, 2, function(x) mean(pdata$ACF*x)/mean(pdata$ACF*!x))
  #Estimates of Productivity Differentials from DS
  DSPB <- apply(XB, 2, function(x) mean(pdata$DS*x)/mean(pdata$DS*!x))
  #Estimates of Productivity Marginal Effects from ACF
  ACFTFPhat <- pdata$ACF[XCind]
  DSTFPhat <- pdata$DS[XCind]
  ACFPC <- apply(XC, 2, function(x) as.numeric(lm(ACFTFPhat~log(x))$coef)[-1])
  #Estimates of Productivity Marginal Effects from DS
  DSPC <- apply(XC, 2, function(x) as.numeric(lm(DSTFPhat~log(x))$coef)[-1])
  return(list(betahat=betahat, ACFhat=ACFhat, qrhat=qrhat, qdifhat=qdifhat, ACFPB=ACFPB, DSPB=DSPB, ACFPC=ACFPC, DSPC=DSPC))
} 
############################################################################################
#Functions for Estimating ACF Coefficients
###########################################################################################
#Function that defines the residuals
ACF_Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi){
  b <- as.matrix(as.numeric(b))
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
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



