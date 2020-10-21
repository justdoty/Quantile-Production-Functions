#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source('gmmq_aux.R')
#Required for QGMM estimation variations
require(GenSA)
require(pracma)
require(dplyr)
###################################################################################
###################################################################################
#This function initializes the estimation procedure which calls finalQLP
###################################################################################
###################################################################################
LP <- function(idvar, timevar, Y, K, L, proxy, b.init=NULL, R=20){
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
  trueboot <- finalQLP(ind=TRUE, data=data, b.init=b.init, gbar=0)
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
    boot.betas[i,] <- finalQLP(ind=boot.indices[[i]], data=data, b.init=b.init, gbar=gbar)$beta
  }
  #Only return the boostrapped estimates and the estimates at the true data
  return(list(boot.betas, beta))
}

###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(ind, data, b.init, gbar){
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
  firststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  phi0 <- as.numeric(coef(first.stage)[1])
  LP_Labor <-  as.numeric(coef(first.stage)[2])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor-phi0
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  #Instruments
  Z <- as.matrix(newdata$Kcon)
  #Use W.init as an initial weighting matrix to get consistent estimates of theta
  Winit <- solve(tau*(1-tau)*crossprod(Z)/nrow(Z))
  #If not specified, starting point is the first stage estimates
  if (is.null(b.init)){
    b.init <- as.numeric(coef(first.stage)[3])
  } 
  #Solve using the above weighting matrix
  if (ncol(Z)>ncol(newdata$Kcon)){
    kinit <- GenSA(par=b.init, fn=LPgo, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=Winit, gbartrue=gbar, lower=0, upper=1, control=list(max.time=5))$par
    W <- solve(var(gbar(theta=kinit, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon,
      fitlagphi=newdata$philag)))
    khat <- GenSA(par=kinit, fn=LPgo, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=W, gbartrue=gbar, lower=0, upper=1, control=list(max.time=5))$par
    gbartrue <- gbar(theta=khat, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fithpi=newdata$phicon,
      lagphi=newdata$philag)
    return(list(beta=c(khat, LP_Labor), gbartrue=gbartrue))
    
  } else if (ncol(Z)==ncol(newdata$Kcon)){
    kinit <- GenSA(par=b.init, fn=LPgo, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=Winit, h=h, gbartrue=gbar, lower=0, upper=1, control=list(max.time=5))$par
    W <- solve(var(gbar(theta=kinit, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon,
      fitlagphi=newdata$philag)))
    khat <- GenSA(par=kinit, fn=LPgo, mX=newdata$Kcon, mlX=newdata$Klag, mZ=Z, fitphi=newdata$phicon, 
      fitlagphi=newdata$philag, mW=W, gbartrue=gbar, lower=0, upper=1, control=list(max.time=5))$par
    gbartrue <- 0
    return(list(beta=c(khat, LP_Labor), gbar=gbartrue))
  }
}
############################################################################################
#This function calculates the residuals used for the moment equations and objective function
###########################################################################################
Lambda <- function(theta, mX, mlX, fitphi, fitlagphi){
  theta <- as.matrix(as.numeric(theta))
  Omega <- fitphi-mX%*%theta[1:ncol(mX)]
  Omega_lag <- fitlagphi-mlX%*%theta[1:ncol(mX)]
  step1 <- lm(Omega~Omega_lag+Omega_lag^2+Omega_lag^3)
  xifit <- resid(step1)
  return(xifit)
} 
##################################################################################
#This function calculates the smoothed sample moment for a given theta
##################################################################################
gbar <- function(theta, mX, mlX, mZ, fitphi, fitlagphi){
  xifit <- Lambda(theta=theta, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi)
  gbar <- as.matrix(colMeans(mZ*repmat(-xifit, 1, ncol(mZ))))
  return(gbar)
}
################################################################################
#QLP Objective Function
################################################################################
LPgo <- function(theta, mX, mlX, fitphi, fitlagphi, mW, gbartrue){
  gbar <- gbar(theta, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi)-gbartrue
  go <- nrow(mZ)*t(g.bar)%*%mW%*%g.bar
  return(go)
} 
##################################################################################
##################################################################################
##################################################################################














