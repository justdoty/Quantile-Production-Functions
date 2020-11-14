# setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('PFQR/FUN/gmmq_aux.R')
source('PFQR/FUN/gmmq.R')
source('PFQR/FUN/1step.R')
source('PFQR/FUN/2step.R')
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
QLP <- function(tau, idvar, timevar, Y, K, L, proxy, h, b.init=NULL, R=20){
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
  if (h==0){
    trueboot <- finalQLP(tau=tau, h=0, ind=TRUE, data=data, b.init=b.init, gbar=0)
  } else {
  	trueboot <- finalQLP(tau=tau, h=h, ind=TRUE, data=data, b.init=b.init, gbar=0)
  }
  #True parameters
  beta <- trueboot$beta
  #True sample moments
  gbar <- trueboot$gbartrue
  #Bandwidth selected from true data
  hbar <- trueboot$h
  #Initialize bootstrap
  boot.indices <- block.boot.resample(idvar, R)
  boot.betas <- matrix(0, nrow=R, ncol=2)
  #Bootstrap Procedure: finalQLP now computes the beta estimates where the sample moments
  #from GMM are recentered by truegbar, the sample moments evaluated at the true data
  for (i in 1:R){
    print(i)
    set.seed(seed+i)
    boot.betas[i,] <- finalQLP(tau=tau, h=hbar, ind=boot.indices[[i]], data=data, b.init=b.init, gbar=gbar)$beta
  }
  #Only return the boostrapped estimates and the estimates at the true data
  return(list(boot.betas, beta, hbar))
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
  phi0 <- as.numeric(coef(firststage)[1])
  LP_Labor <-  as.numeric(coef(firststage)[2])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  print("Data (lag and con)")
  print(head(newdata))
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
  if (is.null(b.init)){
    b.init <- as.numeric(coef(firststage)[3])
  }
  #Overidentification
  if (ncol(mZ)>ncol(mX)){
  	#One Step Estimator of Newey and McFadden (1994)
  	soln <- onestep(mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau, b.init=b.init, gbartrue=gbar, VERBOSE=TRUE)
    #Optional Two-Step Estimator of de Castro et. al. (2019) (uses results from onestep)
    # soln <- twostep(mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, W=soln$W h=soln$h, tau=tau, b.init=soln$b, gbartrue=gbar, upper=1, lower=0, maxtime=5)
    return(list(beta=c(soln$b, LP_Labor), hbar=soln$h, gbartrue=soln$G))
  } else if (ncol(mZ)==ncol(mX)){
    soln <- gmmq(mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau, b.init=b.init, VERBOSE=TRUE)
    gbartrue <- 0
    return(list(beta=c(soln$b, LP_Labor), hbar=soln$h, gbartrue=gbartrue))
  }
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(theta, mY, mX, mlX, fitphi, fitlagphi, tau){
  theta <- as.matrix(as.numeric(theta))
  A <- mY-mX%*%theta[1:ncol(mX)]
  B <- fitlagphi-mX%*%theta[1:ncol(mX)]
  step1 <- rq(A~B+I(B^2)+I(B^3), tau=0.5)
  rho <- as.numeric(coef(step1)[-1])
  residconc <- resid(step1)
  Finv <- quantile(residconc, tau)
  conc <- c(Finv, rho)
  xifit <- A-cbind(1, B, B^2, B^3)%*%conc
  return(xifit)
} 
############################################################################################
#This function calculates the jacobian of Lambda
###########################################################################################
Lambda.derivative <- function(theta, mY, mX, mlX, fitphi, fitlagphi, tau){
  theta <- as.matrix(as.numeric(theta))
  A <- mY-mX%*%theta[1:ncol(mX)]
  B <- fitlagphi-mX%*%theta[1:ncol(mX)]
  step1 <- rq(A~B+I(B^2)+I(B^3), tau=0.5)
  rho <- as.numeric(coef(step1)[-1])
  deriv <- -mX+mlX*rowSums(sweep(cbind(1,B,B^2), margin=2, seq(1,length(rho))*rho, `*`))
  return(deriv)
} 
##################################################################################
##################################################################################
##################################################################################




