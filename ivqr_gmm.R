setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')

#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/

#Some data preparation follows prodest.R (Gabrielle Rovigatti)

#This version implements a Murphy-Topel Estimator for the quantile production function
source("gmmq.R")
require(quantreg)
require(GenSA)
require(pracma)
require(snow)
QACF <- function(tau, va, state, free, proxy, id, time, h=0, b.init=0) {
  #Make all data arguments into matrices
  va <- as.matrix(va)
  state <- as.matrix(state)
  free <- as.matrix(free)
  proxy <- as.matrix(proxy)
  id <- as.matrix(id)
  time <- as.matrix(time)
  #function to generated lagged variables in a panel, see prodest.R (Gabriele Rovigatti) 
  lagPanel <- function(id, time, data){
    df <- data.frame(id, time, data)
    last.time <- df %>% filter(!is.na(time)) %>%
      mutate(time=time + 1, lagged_value=data, data=NULL)
    out <- as.matrix(df %>% left_join(last.time, by = c("id", "time")))[,4]
    return(out)
  }
  #Generate lag values of state and free variable
  lagstate <- state
  for (lag in 1:ncol(state)){
    lagstate[,lag] <- lagPanel(id=id, time=time, data=state[,lag])
  }
  lagfree <- free
  for (lag in 1:ncol(free)){
    lagfree[,lag] <- lagPanel(id=id, time=time, data=free[,lag])
  }
  #Variables to use in first stage regression
  polyframe <- data.frame(state,free,proxy) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(mod, state^2, free^2, proxy^2)

  tmp.data <- as.matrix(data.frame(id=id, time=time, va=va, regvars=regvars))

  firststage <- rq(tmp.data[,'va', drop = FALSE] ~ tmp.data[, grepl('regvars', colnames(tmp.data)), drop = FALSE], na.action = na.exclude, tau=tau)
  gammafirst <- c(as.numeric(coef(firststage)))
  numgamma <- length(gammafirst)
  phi <- fitted(firststage)
  lagphi <- lagPanel(id = id, time=time, data=phi)
  data <- as.matrix(na.omit(data.frame(id=id, time=time, va=data.frame(va), Z=data.frame(lagphi, state, lagfree), 
      Xt=data.frame(state, free), lX=data.frame(lagstate, lagfree), Lagphi=data.frame(lagphi), regvars=regvars)))
  va <- data[,"va"]
  Z <- data[, grepl('Z', colnames(data)), drop = FALSE]
  Xt <- data[, grepl('Xt', colnames(data)), drop = FALSE]
  lX <- data[, grepl('lX', colnames(data)), drop = FALSE]
  Lagphi <- data[,"lagphi"]
  regvars <- data[, grepl('regvars', colnames(data)), drop = FALSE]

  #Sample size
  n <- length(id)
  #Number of instruments
  dZ <- ncol(Z)
  Lambda <- function(va, Xt, lX, Lagphi, b){
    Moment <- va-Xt%*%b[1:(length(b)-1)]-b[length(b)]*(Lagphi-lX%*%b[1:(length(b)-1)])
    return(Moment)
  }
  #Smoothing Kernels
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
  # #Define Moment Equation

  # #Add bandwidth as second argument to G() and G'() functions.
  Gfn <- function(v,h){      
    Itilde.KS17(v/h)    
  }
  Gpfn <- function(v,h){      
    Itilde.deriv.KS17(v/h)    
  }
######### Second Step: Initial Consistent Estimates and Mu Estimates ########################

  # #Iniitial arbitrary weight matrix
  W.init <- diag(dZ)
  # #Second step estimation objective function
  ivqr.obj <- function(b, h, W) {
      L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b), ncol=1)
      gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
      g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 
      obj.fn <- t(g.bar)%*%W%*%g.bar
      return(obj.fn)
  }
  obj.fn <- function(b, W) {
      return(n*ivqr.obj(b, h=h, W))
  }
  ##Second Step Moments
  M.hat <- function(b) {
    L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b), ncol=1)
    gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
    g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 
    return(g.bar)        
  }
  ##Minimize the sample moment with the initial arbitrary weight matrix W.init 
  ivqr.gmm1 <- GenSA(par=b.init, fn=function(b) {obj.fn(b, W=W.init)}, lower=array(0, dZ), upper=array(1, dZ), control=list(max.time=20))
  # #Initial consistent estimator
  b1 <- ivqr.gmm1$par
  ########## Second Step: Initial Consistent Estimates and Mu Estimates ########################
  ########## ACHL Second Step ##########################################
  ########## Get estimates of the mu's
  mu.hat <- as.numeric(M.hat(b1))
  ########## Collect New Parameters
  gammamu <- c(gammafirst, mu.hat)
  #Joint Moments
  joint.mom <- function(locgammamu){
    lhhat <- cbind(1, regvars)%*%locgammamu[1:length(gammafirst)]
    resid1 <- va-lhhat
    resid2 <- va-Xt%*%b1[1:(length(b1)-1)]-b1[length(b1)]*(lhhat-lX%*%b1[1:(length(b1)-1)])
    mom1 <- cbind(1,regvars)*repmat((Gfn(-resid1,h)-tau),1,(ncol(regvars)+1))
    mom2 <- Z*repmat((Gfn(-resid2,h)-tau),1,dZ)-locgammamu[(length(gammafirst)+1):length(gammamu)] 
    mom <- cbind(mom1, mom2)
    return(mom)
  }
  # Derivative of Joint Moments
  joint.mom.der <- function(locgammamu){
    lhhat <- cbind(1, regvars)%*%locgammamu[1:length(gammafirst)]
    L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=lhhat, b=b1), ncol=1)
    L1 <- va-lhhat
    L2 <- va-Xt%*%b1[1:(length(b1)-1)]-b1[length(b1)]*(lhhat-lX%*%b1[1:(length(b1)-1)])
    L1.Step1 <- array(data=Itilde.deriv(-L1/h),dim=dim(cbind(1, regvars)))*cbind(1,regvars)
    L2.Step1 <- array(data=Itilde.deriv(-L2/h),dim=dim(Z))*Z
    L1.d <- -cbind(1, regvars)
    L2.d <- -locgammamu[length(gammafirst)]*cbind(1, regvars)
    G11 <- -t(L1.Step1)%*%L1.d
    G21 <- -t(L2.Step1)%*%L2.d
    #Construct the Entire G Matrix by Blocks
    G12 <- matrix(0, nrow=(ncol(regvars)+1), ncol=length(mu.hat))
    G22 <- matrix(-1, nrow=length(mu.hat), ncol=length(mu.hat))
    G <- cbind(rbind(G22, G12), rbind(G21, G11))
    return(G)
  }
  #Now do Murphy-Topel Estimate of Variance of Mu
  varmat <- solve(cov(joint.mom(gammamu)))
  G <- joint.mom.der(gammamu)
  avar <- pinv(G%*%varmat%*%t(G))*n
  weight.mat <- avar[1:length(mu.hat), 1:length(mu.hat)]
  ivqr.gmm2 <- GenSA(par=b1, fn=function(b) {obj.fn(b, W=weight.mat)}, lower=array(0, dZ), upper=array(1, dZ), control=list(max.time=5))
  b2 <- ivqr.gmm2$par
  return(b1)
}
hi <- QACF(tau=0.9, va=chile$lnva, state=chile$lnk, free=cbind(chile$lnb, chile$lnw), proxy=chile$proxy_m, id=chile$ppn, time=chile$year, h=0.001, b.init=c(0,0,0,0))



