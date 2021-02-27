# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('PFQR/FUN/QLP_aux.R')
source('PFQR/FUN/BW.R')
require(quantreg)
require(dplyr)
require(pracma)
require(GenSA)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
QLP <- function(tau, h, idvar, timevar, Y, K, L, proxy, dZ, binit=NULL){
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
  res <- finalQLP(tau=tau, h=h, data=data, binit=binit, seed=seed)
  #"True" parameters
  betahat <- res$beta
  LPhat <- res$LPhat
  return(list(betahat=betahat, LPhat=LPhat))
}
###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(tau, h, data, binit, seed){
  set.seed(seed)
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
  wfit <- newdata$LPphicon-newdata$Kcon*LPkhat
  #Output net of productivity
  mY <- as.matrix(newdata$Ycon-wfit)
  #State Variables
  mX <- cbind(newdata$Kcon, newdata$Lcon)
  #Instruments
  mZ <- cbind(as.matrix(newdata$Kcon), as.matrix(newdata$Lcon), as.matrix(newdata$Klag), as.matrix(newdata$Llag), as.matrix(newdata$Pxcon), as.matrix(newdata$Pxlag))
  #QLP Estimates for Capital from (QR) good starting values as well
  mom <- rq(mY~mX-1, tau=tau)
  init <- as.numeric(coef(mom))
  if (is.null(h)){
    h <- ivqr.bw(p=tau, Y=mY, X=mX, b.init=init)
  }
  betahat <- GenSA(par=init, fn=QLPobj, mY=mY, mX=mX, mZ=mZ, h=h, tau=tau, lower=c(0,0), upper=c(1,1), control=list(max.time=1))$par
  #LP Estiates
  LPhat <- c(LPkhat, LPLabor)
  return(list(betahat=betahat, LPhat=LPhat))
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(b, mY, mX, Qfit){
  resid <- mY-mX*b-Qfit
  return(resid)
} 
#QLP GMM objective function
 QLPobj <- function(b, mY, mX, mZ, Qfit, tau, h){
  resid <- Lambda(b=b, mY=mY, mX=mX, Qfit=Qfit)
  gni <- mZ*repmat((Gfn(-resid, h)-tau), 1, ncol(mZ))
  gbar <- as.matrix(colMeans(gni))
  W <- solve(tau*(1-tau)*t(mZ)%*%mZ/nrow(mZ))
  go <- nrow(mZ)*t(gbar)%*%W%*%gbar
  return(go)

} 
############################################################################################
#Functions for Estimating QLP Coefficients
###########################################################################################
#Function that defines the residuals
Lambda <- function(b, mY, mX){
  b <- as.matrix(as.numeric(b))
  resid <- mY-mX%*%b
  return(resid)
} 
#QLP GMM objective function
 QLPobj <- function(b, mY, mX, mZ, gbar, tau, h){
  resid <- Lambda(b=b, mY=mY, mX=mX)
  gni <- mZ*repmat((Gfn(-resid, h)-tau), 1, ncol(mZ))
  gbar <- as.matrix(colMeans(gni))
  W <- solve(tau*(1-tau)*t(mZ)%*%mZ/nrow(mZ))
  go <- nrow(mZ)*t(gbar)%*%W%*%gbar
  return(go)
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
LPobj <- function(b, mY, mX, mZ, mlX, fitphi, fitlagphi){
  resid <- LP_Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)
  gni <- mZ*repmat(resid,1, ncol(mZ))
  gnic <- colSums(gni)
  go <- sum(gnic^2)
  return(go)
}
