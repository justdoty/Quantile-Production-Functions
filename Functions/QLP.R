




































     gbar <- as.matrix(colMeans(gni)) 
     gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
     return(gbar)
     xifit <- matrix(Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
    binit <- as.numeric(coef(firststage)[3])
  #####################
  #Calculate Contempory and Lag Values for 2nd stage estimation
  #Clean Phi from the effects of free variables
  #Contemporary phi estimates
  #Contemporary State Variables
  #Instruments
  #Lagged phi estimates
  #Lagged State Variables
  #Make all data arguments into matrices
  #Output net of labor
  5
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  betal <-  firststage$coefficients[2,1]
  betalCI <- firststage$coefficients[2,2:3]
  betalse <- summary(firststage)$coefficients[2,2]
  data <- data.frame(idvar=idvar, timevar=timevar, Y=Y, K=K, L=L, proxy=proxy)
  firststage <- rq(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]), tau=tau, alpha=alpha, ci=TRUE)
  fitlagphi <- as.matrix(newdata$philag)
  fitphi <- as.matrix(newdata$phicon)
  gbar <- as.matrix(colMeans(gni))
  gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
  go <- nrow(mY)*t(gbar)%*%W%*%gbar
  idvar <- as.matrix(idvar)
  if (is.null(binit)){
  K <- as.matrix(K)
  L <- as.matrix(L)
  mlX <- as.matrix(newdata$Klag)
  mX <- as.matrix(newdata$Kcon)
  mY <- as.matrix(newdata$Ycon-newdata$Lcon*LP_Labor)
  mZ <- cbind(as.matrix(newdata$Kcon), as.matrix(newdata$Klag), as.matrix(newdata$Llag))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor
  proxy <- as.matrix(proxy)
  regvars <- data.frame(reg1=data$L, reg2=data$K, reg3=data$proxy, reg4=data$K^2, reg5=data$proxy^2)
  return(go)
  return(list(boot.betas, beta))
  return(xifit)
  seed <- 123456
  set.seed(seed)
  step1 <- lm(A~B+I(B^2)+I(B^3))
  step1param <- as.numeric(coef(step1))
  timevar <- as.matrix(timevar)
  W <- solve(var(gni))
  xifit <- A-cbind(1, B, B^2, B^3)%*%step1param
#   deriv <- -mX+mlX*rowSums(sweep(cbind(1,B,B^2), MARGIN=2, seq(1,length(rho))*rho, `*`))
#   return(deriv)
#   rho <- as.numeric(coef(step1))
#   step1 <- lm(A~B+I(B^2)+I(B^3))
# Lambda.derivative <- function(b, mY, mX, mlX, fitphi, fitlagphi, tau){
# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
# } 
##################################################################################
##################################################################################
##################################################################################
###################################################################################
###################################################################################
###########################################################################################
###########################################################################################
############################################################################################
############################################################################################
#https://faculty.missouri.edu/~kaplandm/
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#Required for 1st step estimation
#Required for QGMM estimation variations
#Some data preparation follows prodest.R (Gabrielle Rovigatti)
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#This function calculates the jacobian of Lambda
#This function calculates the residuals used for the moment equations and objective function
gbar <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi, gbar, tau){
Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi, tau){
QLP <- function(tau, idvar, timevar, Y, K, L, proxy, binit=NULL, alpha){
require(dplyr)
require(pracma)
require(quantreg)
source('PFQR/FUN/gmmq_aux.R')
}
}
} 