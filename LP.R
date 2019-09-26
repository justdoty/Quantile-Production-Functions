setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
set.seed(123456)
#This version implements a Murphy-Topel Estimator for Ackerberg Caves Frazer (2015)
#"Identification Properties of Recent Production Function Estimators" as in Ackerberg Chen, Hahn (2012)"
#"A Practical Asymptotic Variance Estimator for Two-Step Semiparametric Estimators"

#Some data preparation follows prodest.R (Gabrielle Rovigatti)
require(GenSA)
require(dplyr)
require(pracma)
require(MASS)
require(prodest)
LP <- function(va, state, free, proxy, id, time, b.init=NULL){
  # Make all data arguments into matrices
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
  polyframe <- data.frame(state, proxy) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(free, mod, state^2, proxy^2)

  tmp.data <- suppressWarnings(as.matrix(data.frame(id=id, time=time, va=va, regvars=regvars)))
  firststage <- lm(tmp.data[,'va', drop = FALSE] ~ tmp.data[, grepl('regvars', colnames(tmp.data)), drop = FALSE], na.action = na.exclude)
  gammafirst <- c(as.numeric(firststage$coefficients))
  gammal <- gammafirst[2:(ncol(free)+1)]
  numgamma <- length(gammafirst)
  phi <- fitted(firststage)-regvars[,1:ncol(free)]%*%as.matrix(gammal)
  lagphi <- lagPanel(id = id, time=time, data=phi)
  data <- suppressWarnings(as.matrix(na.omit(data.frame(id=id, time=time, va=data.frame(va), Z=data.frame(lagstate, lagphi), 
        state=data.frame(state), lagstate=data.frame(lagstate), lagphi=data.frame(lagphi), regvars=regvars))))
  va <- as.matrix(data[,"va"])
  Z <- data[, grepl('Z', colnames(data)), drop = FALSE]
  Xt <- as.matrix(data[,"state"])
  lX <- as.matrix(data[, "lagstate"])
  Lagphi <- as.matrix(data[,"lagphi"])
  regvars <- data[, grepl('regvars', colnames(data)), drop = FALSE]
  va.net <- va-regvars[,1:ncol(free)]%*%as.matrix(gammal)
  #Sample size
  n <- length(va)
  ##Number of instruments
  dZ <- ncol(Z)
  LP_GMM <- function(va, Xt, lX, Lagphi, Z, b){
	  Moment <- va-Xt%*%b[1:ncol(Xt)]-b[length(b)]*(Lagphi-lX%*%b[1:ncol(lX)])
	  Obj <- Z*array(data=Moment, dim=dim(Z))
	  return(Obj)
  } 
  obj.fn <- function(b, W){
  	momi <- LP_GMM(va=va.net, Xt=Xt, lX=lX, Lagphi=Lagphi, Z=Z, b)
  	return(n*colMeans(momi)%*%W%*%as.matrix(colMeans(momi)))
  }
  M.hat <- function(b){
  	momi <- LP_GMM(va=va.net, Xt=Xt, lX=lX, Lagphi=Lagphi, Z=Z, b)
  	return(colMeans(momi))
  }
  W.init <- diag(dZ)
  # # # stage1 <- GenSA(par=b.init, fn=function(b){obj.fn(b, W=W.init)}, lower=array(0, dZ), upper=array(1, dZ), control=list(max.time=5))
  if (is.null(b.init)){
  	b.init <- c(gammafirst[(ncol(free)+2):(ncol(free)+ncol(state)+1)], 0.7)+rnorm(dZ, 0, 0.01)
  } 
  stage1 <- suppressWarnings(solnp(pars=b.init, fun=function(b){obj.fn(b, W=W.init)}, control=list(trace=FALSE)))
  b1 <- stage1$pars
  mu.hat <- as.numeric(M.hat(b1))
  gammamu <- c(gammafirst, mu.hat)

  joint.mom <- function(locgammamu){
  	#Fitted Labor
  	lhhat <- regvars[,1:ncol(free)]%*%as.matrix(locgammamu[2:(ncol(free)+1)])
  	#Fitted Phi net of fitted labor
    phihat <- cbind(1, regvars)%*%locgammamu[1:length(gammafirst)]-lhhat
    resid1 <- va-lhhat-phihat
    mom1 <- cbind(1,regvars)*repmat(resid1, 1,(ncol(regvars)+1))
    mom2 <- LP_GMM(va=va.net, Xt=Xt, lX=lX, Lagphi=phihat, Z=Z, b=b1)-locgammamu[(length(gammafirst)+1):length(gammamu)] 
    mom <- cbind(mom1, mom2)
    return(mom)
  }
  ###Derivative of Joint Moments
  # joint.mom.der <- function(locgammamu){
  # 	G11 <- -diag(length(mu.hat))
  # 	G21 <- matrix(0, nrow=(ncol(regvars)+1), ncol=length(mu.hat))
  # 	G12.Free <- -t(Z)%*%regvars[,1:ncol(free)]
  # 	G12.Phi <- -locgammamu[length(mu.hat)]*t(Z)%*%cbind(1, regvars[,-(1:ncol(free))])
  # 	G12 <- cbind(G12.Free, G12.Phi)/n
  # 	G22.Free <- -t(cbind(1, regvars))%*%regvars[,1:ncol(free)]
  # 	G22.Phi <-  -t(cbind(1, regvars))%*%cbind(1, regvars[,-(1:ncol(free))])
  # 	G22 <- cbind(G22.Free, G22.Phi)/n
  # 	G <- rbind(cbind(G22, G21), cbind(G12, G11))
  # 	return(G)
  # }
  #Optional for numeric derivatives
  joint.deriv <- function(locgammamu){
  	return(colMeans(joint.mom(locgammamu)))
  }
  #Numeric jacobian seems to perform better than analytical
  G.numeric <- jacobian(joint.deriv, gammamu)
  # # # Now do Murphy-Topel Estimate of Variance of Mu
  # # # The Estimate of the covariance of the joint moment restrictions
  varmat <- crossprod(joint.mom(gammamu))/n
  # # # The estimate of the derivatives of the joint moment restrictions
  G <- solve(G.numeric)
  # G <- solve(joint.mom.der(gammamu))
  # # # #The estimate of the asymptotic covariance matrix of gammas and mu's
  avar <- G%*%varmat%*%t(G)/n
  # # # #Take the block corresponding to the covariance matrix of mu
  weight.mat <- solve(avar[(nrow(avar)-dZ+1):nrow(avar), (ncol(avar)-dZ+1):ncol(avar)])
  #Take the block corresponding to the covariance matrix of free variables
  if (!is.matrix(avar[1:ncol(free), 1:ncol(free)])){
      avar.l <- as.numeric(diag(avar[1:ncol(free), 1:ncol(free)], nrow=1, ncol=1))
    } else {
      avar.l <- as.numeric(diag(avar[1:ncol(free), 1:ncol(free)]))
    }
  # # # # #Take the inverse for the weighting matrix to use in the last step
  stage2 <- suppressWarnings(solnp(pars=b1, fun=function(b){obj.fn(b, W=weight.mat)}, control=list(trace=FALSE)))
  b2 <- stage2$pars
  # # #Now calculate the asymptotic variance of ACF_GMM.
  Moment.derivative <- function(b){
  	#State Variable
    g1 <- -Xt[,1:ncol(state)]+b[length(b)]*lX[,1:ncol(state)]
    # #Productivity
    g2 <- -(Lagphi-lX%*%b[1:(length(b)-1)])
    ##Combined
    g <- t(Z)%*%cbind(g1, g2, deparse.level = 0)/n
    return(g)
  }
  G.b <- solve(Moment.derivative(b2))
  avar.b <- as.numeric(diag(G.b%*%weight.mat%*%t(G.b)/n))
  return(list(c(b2[1:ncol(state)], gammal, b2[length(b2)]), c(avar.b[1:ncol(state)], avar.l, avar.b[length(avar.b)])))
}

# chile_panel <- read.csv('chile_panel.csv')
# chile <- na.omit(subset(chile_panel, ciiu_3d==381))

# results <- LP(va=chile$lnva, state=chile$lnk, free=chile$lnl, proxy=chile$proxy_e, id=chile$id, time=chile$year)
# print(results)