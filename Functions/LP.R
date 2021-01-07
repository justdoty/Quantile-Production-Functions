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
LP <- function(idvar, timevar, Y, K, L, proxy, binit=NULL, R=20){
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
  trueboot <- finalQLP(ind=TRUE, data=data, binit=binit, gbar=0, seed=seed)
  #"True" parameters
  betahat <- trueboot$beta
  #"True" TFP dispersion rations
  ratiohat <- trueboot$dispersion
  #True sample moments
  gbartrue <- trueboot$gbar
  #Initialize bootstrap
  bootind <- block.boot.resample(idvar, R, seed)
  betaboot <- matrix(0, nrow=R, ncol=2)
  ratioboot <- matrix(0, nrow=R, ncol=3)
  #Bootstrap Procedure: finalQLP now computes the beta estimates where the sample moments
  #from GMM are recentered by truegbar, the sample moments evaluated at the true data
  for (i in 1:R){
    print(i)
    seed <- seed+i
    boot <- finalQLP(ind=bootind[[i]], data=data, binit=binit, gbar=gbartrue, seed=seed)
    betaboot[i,] <- boot$beta
    ratioboot[i,] <- boot$dispersion
  }
  return(list(betahat=betahat, ratiohat=ratiohat,  betaboot=betaboot, ratioboot=ratioboot))
}

###########################################################################
###########################################################################
#Function to estimate and to bootstrap QLP
###########################################################################
###########################################################################
finalQLP <- function(ind, data, binit, gbar, seed){
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
  #Create a polynomail data frame for 1st stage estimation
  regvars  <- data.frame(reg1=data$L, reg2=data$K, reg3=data$proxy, reg4=data$proxy*data$K, reg5=data$K^2, reg6=data$proxy^2, reg7=data$proxy*(data$K^2), reg8=(data$proxy^2)*data$K, reg9=data$K^3, reg10=data$proxy^3)
  firststage <- lm(data$Y~as.matrix(regvars[, grepl('reg', colnames(regvars))]))
  init <- lm(data$Y~data$L+data$K)
  kinit <- as.numeric(coef(init)[3])
  phi0 <- as.numeric(coef(firststage)[1])
  LP_Labor <-  as.numeric(coef(firststage)[2])
  #Clean Phi from the effects of free variables
  phi <- fitted(firststage)-as.matrix(data$L)%*%LP_Labor
  #Calculate Contempory and Lag Values for 2nd stage estimation
  newdata <- lagdata(idvar=data$idvar, X=cbind(data$Y, data$K, data$L, data$proxy, phi))
  names(newdata) <- c("idvar", "Ycon", "Kcon", "Lcon", "Pxcon", "phicon", "Ylag", "Klag", "Llag", "Pxlag", "philag")
  #Output net of labor
  mY <- as.matrix(newdata$Ycon-newdata$Lcon*LP_Labor)
  mX <- as.matrix(newdata$Kcon)
  mlX <- as.matrix(newdata$Klag)
  fitphi <- as.matrix(newdata$phicon)
  fitlagphi <- as.matrix(newdata$philag)
  #Instruments
  mZ <- cbind(as.matrix(newdata$Kcon), as.matrix(newdata$Klag), as.matrix(newdata$Llag), as.matrix(newdata$Pxlag))
  #If not specified, starting point is the first stage estimates
  if (is.null(binit)){
     # binit <- 1-LP_Labor
     binit <- kinit
  } 
  if (ncol(mZ)>ncol(mX)){
    soln <- optim(par=binit, fn=function(b) LPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    gbar <- gbar(b=soln$par, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar)
    TFP <- exp(data$Y-cbind(data$K, data$L)%*%c(soln$par, LP_Labor))
    Q3Q1hat <- as.numeric(quantile(TFP, .75)/quantile(TFP, .25))
    Q9Q1hat <- as.numeric(quantile(TFP, .9)/quantile(TFP, .1))
    Q95Q05hat <- as.numeric(quantile(TFP, .95)/quantile(TFP, .05))
    betahat <- c(soln$par, LP_Labor)
    dispersion <- c(Q3Q1hat, Q9Q1hat, Q95Q05hat)
    print(head(data))
    print(betahat)
    print(dispersion)
    return(list(betahat=betahat, gbar=gbar, dispersion=dispersion))
    
  } else if (ncol(mZ)==ncol(mX)){
    soln <- optim(par=binit, fn=function(b) LPobj(b, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, gbar=gbar), gr=NULL, method="L-BFGS-B", lower=0, upper=1)
    gbar <- 0
    TFP <- exp(data$Y-cbind(data$K, data$L)%*%c(soln$par, LP_Labor))
    Q3Q1hat <- as.numeric(quantile(TFP, .75)/quantile(TFP, .25))
    Q9Q1hat <- as.numeric(quantile(TFP, .9)/quantile(TFP, .1))
    Q95Q05hat <- as.numeric(quantile(TFP, .95)/quantile(TFP, .05))
    betahat <- c(soln$par, LP_Labor)
    dispersion <- c(Q3Q1hat, Q9Q1hat, Q95Q05hat)
    print(betahat)
    print(dispersion)
    return(list(betahat=betahat, gbar=gbar, dispersion=dispersion))
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














