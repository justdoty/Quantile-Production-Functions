# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('PFQR/FUN/gmmq_aux.R')
#Required for 1st step estimation
require(quantreg)
#Required for QGMM estimation variations
require(pracma)
require(dplyr)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
QLP <- function(tau, idvar, timevar, Y, K, L, proxy, binit=NULL, R=20){
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
  trueboot <- finalQLP(tau=tau, ind=TRUE, data=data, binit=binit, gbar=0)
  #True parameters
  beta <- trueboot$beta
  #True sample moments
  gbartrue <- trueboot$gbar
  #Initialize bootstrap
  boot.indices <- block.boot.resample(idvar, R)
  boot.betas <- matrix(0, nrow=R, ncol=2)
  #Bootstrap Procedure: finalQLP now computes the beta estimates where the sample moments
  #from GMM are recentered by truegbar, the sample moments evaluated at the true data
  for (i in 1:R){
    print(i)
    set.seed(seed+i)
    boot.betas[i,] <- finalQLP(tau=tau, ind=boot.indices[[i]], data=data, binit=binit, gbar=gbartrue)$beta
  }
  #Only return the boostrapped estimates and the estimates at the true data
  return(list(boot.betas, beta))
}
###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(tau, ind, data, binit, gbar){
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
  regvars <- data.frame(reg1=data$L, reg2=data$K, reg3=data$proxy, reg4=data$K^2, reg5=data$proxy^2)
  firststage <- rq(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]), tau=tau)
  phi0 <- as.numeric(coef(firststage)[1])
  LP_Labor <-  as.numeric(coef(firststage)[2])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  #Output net of labor
  mY <- as.matrix(newdata$Ycon-newdata$Lcon*LP_Labor)
  #Contemporary State Variables
  mX <- as.matrix(newdata$Kcon)
  #Lagged State Variables
  mlX <- as.matrix(newdata$Klag)
  #Contemporary phi estimates
  fitphi <- as.matrix(newdata$phicon)
  #Lagged phi estimates
  fitlagphi <- as.matrix(newdata$philag)
  #Instruments
  mZ <- cbind(as.matrix(newdata$Kcon), as.matrix(newdata$Klag), as.matrix(newdata$Llag))
  #Optional weighting matrix used by Kaplan and Sun (2016)
  # Wz <- solve(tau*(1-tau)*crossprod(mZ)/nrow(mZ))
  #If not specified, starting point is the first stage estimates
  if (is.null(binit)){
    binit <- as.numeric(coef(firststage)[3])
  }
  #Overidentification
  if (ncol(mZ)>ncol(mX)){
  	# soln <- GenSA(par=binit, fn=QLPobj, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar, tau=tau, lower=0, upper=1, control=list(max.time=5))
    soln <- optim(par=binit, fn=function(b) QLPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar, tau=tau), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    gbar <- gbar(b=soln$par, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar, tau=tau)
    return(list(beta=c(soln$par, LP_Labor), gbar=gbar))
  } else if (ncol(mZ)==ncol(mX)){
    soln <- optim(par=binit, fn=function(b) QLPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar, tau=tau), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    gbar <- 0
    return(list(beta=c(soln$b, LP_Labor), gbar=gbar))
  }
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(b, mY, mX, mlX, fitphi, fitlagphi, tau){
  A <- fitphi-mX%*%b[1:ncol(mX)]
  B <- fitlagphi-mlX%*%b[1:ncol(mX)]
  step1 <- lm(A~B+I(B^2)+I(B^3))
  step1param <- as.numeric(coef(step1))
  xifit <- A-cbind(1, B, B^2, B^3)%*%step1param
  return(xifit)
} 
gbar <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi, gbar, tau){
     xifit <- matrix(Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
     gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
     gbar <- as.matrix(colMeans(gni)) 
     return(gbar)
  }
 QLPobj <- function(b, mY, mX, mlX, mZ, fitphi, fitlagphi, gbar, tau){
  xifit <- Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)
  gni <- sweep(mZ*repmat(xifit,1, ncol(mZ)), MARGIN=2, gbar, `-`)
  gbar <- as.matrix(colMeans(gni))
  W <- solve(var(gni))
  go <- nrow(mY)*t(gbar)%*%W%*%gbar
  return(go)

}
############################################################################################
#This function calculates the jacobian of Lambda
###########################################################################################
# Lambda.derivative <- function(b, mY, mX, mlX, fitphi, fitlagphi, tau){
#   b <- as.matrix(as.numeric(b))
#   A <- fitphi-mX%*%b[1:ncol(mX)]
#   B <- fitlagphi-mlX%*%b[1:ncol(mX)]
#   step1 <- lm(A~B+I(B^2)+I(B^3))
#   rho <- as.numeric(coef(step1))
#   deriv <- -mX+mlX*rowSums(sweep(cbind(1,B,B^2), MARGIN=2, seq(1,length(rho))*rho, `*`))
#   return(deriv)
# } 
##################################################################################
##################################################################################
##################################################################################




