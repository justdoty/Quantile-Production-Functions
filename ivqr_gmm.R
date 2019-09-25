setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
set.seed(123456)
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/

#Some data preparation follows prodest.R (Gabrielle Rovigatti)

#This version implements a Murphy-Topel Estimator for the quantile production function
source("gmmq.R")
source("ivqr.bw.R")
#Required for 1st step estimation
require(quantreg)
#Required for QGMM estimation
require(GenSA)
require(pracma)
#Optional for parallel computing
require(snow)
require(dplyr)
#Required for IVQR.BW
require(MASS)
#Required for ACF estimation comparison
require(prodest)
QACF <- function(tau, va, state, free, proxy, id, time, h=0, b.init=0){
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

  tmp.data <- suppressWarnings(as.matrix(data.frame(id=id, time=time, va=va, regvars=regvars)))

  #Estimate firststage parameters (gamma) and fit model
  firststage <- rq(tmp.data[,'va', drop = FALSE] ~ tmp.data[, grepl('regvars', colnames(tmp.data)), drop = FALSE], na.action = na.exclude, tau=tau)
  gammafirst <- c(as.numeric(coef(firststage)))
  numgamma <- length(gammafirst)
  phi <- fitted(firststage)
  lagphi <- lagPanel(id = id, time=time, data=phi)
  data <- suppressWarnings(as.matrix(na.omit(data.frame(id=id, time=time, va=data.frame(va), Z=data.frame(state, lagfree, lagphi), 
        Xt=data.frame(state, free), lX=data.frame(lagstate, lagfree), Lagphi=data.frame(lagphi), regvars=regvars))))
  va <- data[,"va"]
  Z <- data[, grepl('Z', colnames(data)), drop = FALSE]
  Xt <- data[, grepl('Xt', colnames(data)), drop = FALSE]
  lX <- data[, grepl('lX', colnames(data)), drop = FALSE]
  Lagphi <- data[,"lagphi"]
  regvars <- data[, grepl('regvars', colnames(data)), drop = FALSE]
  ##Sample size
  n <- length(va)
  ##Number of instruments
  dZ <- ncol(Z)
  #Moment function in 2nd step (plugging in fitted firststage values)
  Lambda <- function(va, Xt, lX, Lagphi, b){
    Moment <- va-Xt%*%b[1:(length(b)-1)]-b[length(b)]*(Lagphi-lX%*%b[1:(length(b)-1)])
    return(Moment)
  }
  #Derivative of moment function in 2nd step (plugging in fitted firststage values)
  Lambda.derivative <- function(va, Xt, lX, Lagphi, b){
    #State Variables
    g1 <- -Xt[,1:ncol(state)]+b[length(b)]*lX[,1:ncol(state)]
    #Free Variables
    g2 <- -Xt[,(ncol(state)+1):ncol(Xt)]+b[length(b)]*lX[,(ncol(state)+1):ncol(Xt)]
    # #Productivity
    g3 <- -(Lagphi-lX%*%b[1:(length(b)-1)])
    # # #Combined
    g <- cbind(g1, g2, g3, deparse.level = 0)
    return(g)
    }
    #Add bandwidth as second argument to G() and G'() functions.
    Gfn <- function(v,h){      
    Itilde.KS17(v/h)    
    }
    Gpfn <- function(v,h){      
      Itilde.deriv.KS17(v/h)    
    }
  #Smoothing Kernels
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
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
  ##Second Step Sample Moments
  M.hat <- function(b) {
    L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b), ncol=1)
    gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
    g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 
    return(g.bar)  
  } 
########First Stage: Initial Consistent Estimates#########################
  #Here i use the over identified case for estimates of beta using iid weighting matrix of Kaplan and Sun (2017)
  W.init <- solve((1-tau)*tau*t(Z)%*%Z/n)
  ivqr.est <- GenSA(par=b.init, fn=function(b){obj.fn(b, W=W.init)}, lower=array(0, dZ), upper=array(1, dZ), control=list(max.time=5))
  b1 <- ivqr.est$par
  #Below implements a first stage consistent estimate using gmmq (exact identification) but does not work well
    # ivqr.est <- gmmq(tau=tau, va=va, Xt=Xt, lX=lX, Z=Z, Lagphi=Lagphi, id=id, time=time, 
    #   Lambda=Lambda, Lambda.derivative=Lambda.derivative, h=0, VERBOSE=TRUE, b.init=b.init)
    # #First stage consistent estimates of beta 
    # b1 <- ivqr.est$b
    #First stage bandwidth
    #   h <- ivqr.est$h
########## Second Stage: Mu Estimates ########################      
########### Get estimates of the mu's for fixed value of beta
  mu.hat <- as.numeric(M.hat(b1))
########### Collect New Parameters
  gammamu <- c(gammafirst, mu.hat)
##Joint Moments from 1st and 2nd step
  joint.mom <- function(locgammamu){
    lhhat <- cbind(1, regvars)%*%locgammamu[1:length(gammafirst)]
    resid1 <- va-lhhat
    resid2 <- va-Xt%*%b1[1:(length(b1)-1)]-b1[length(b1)]*(lhhat-lX%*%b1[1:(length(b1)-1)])
    mom1 <- cbind(1,regvars)*repmat((Gfn(-resid1,h)-tau),1,(ncol(regvars)+1))
    mom2 <- Z*repmat((Gfn(-resid2,h)-tau),1,dZ)-locgammamu[(length(gammafirst)+1):length(gammamu)] 
    mom <- cbind(mom1, mom2)
    return(mom)
  }
##Control for smoothing in G function: The smallest feasible bandwidth is too small for estimating G
    #Plug-in bandwidth from Kaplan and Sun (2017)
    # h.joint <- ivqr.bw(p=tau, va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b.init=b1)
  h.joint <- n^(-1/2)
##Derivative of Joint Moments
  joint.mom.der <- function(locgammamu){
    lhhat <- cbind(1, regvars)%*%locgammamu[1:length(gammafirst)]
    L <- matrix(Lambda(va=va, Xt=Xt, lX=lX, Lagphi=lhhat, b=b1), ncol=1)
    L1 <- va-lhhat
    L2 <- va-Xt%*%b1[1:(length(b1)-1)]-b1[length(b1)]*(lhhat-lX%*%b1[1:(length(b1)-1)])
    L1.Step1 <- array(data=Itilde.deriv(-L1/h.joint),dim=dim(cbind(1, regvars)))*cbind(1,regvars)
    L2.Step1 <- array(data=Itilde.deriv(-L2/h.joint),dim=dim(Z))*Z
    #Derivative of the first moment residual with respect to gammas
    L1.d <- -cbind(1, regvars)
    #Derivative of the second moment residual with respect to gammas
    L2.d <- -locgammamu[length(gammafirst)]*cbind(1, regvars)
    #Construct the Entire G Matrix by Blocks
    G11 <- -diag(length(mu.hat))
    G21 <- matrix(0, nrow=(ncol(regvars)+1), ncol=length(mu.hat))
    G12 <- -(t(L2.Step1)%*%L2.d)
    G22 <- -(t(L1.Step1)%*%L1.d)
    G <- cbind(rbind(G11, G21), rbind(G12, G22))/(n*h.joint)
    return(G)
  }
  #Now do Murphy-Topel Estimate of Variance of Mu
  #The Estimate of the covariance of the joint moment restrictions
  varmat <- t(joint.mom(gammamu))%*%joint.mom(gammamu)/n
  # The estimate of the derivatives of the joint moment restrictions
  G <- solve(joint.mom.der(gammamu))
  #The estimate of the asymptotic covariance matrix of gammas and mu's
  avar <- G%*%varmat%*%t(G)
  #Take the block corresponding the the covariance matrix of mu
  var <- avar[1:length(mu.hat), 1:length(mu.hat)]/n
  # #Take the inverse for the weighting matrix to use in the last step
  weight.mat <- solve(var)
  ivqr.gmm2 <- GenSA(par=b1, fn=function(b){obj.fn(b, W=weight.mat)}, lower=array(0, dZ), upper=array(1, dZ), control=list(max.time=5))
  b2 <- ivqr.gmm2$par
  # Derivative of Moments for b
  #Plug-in bandwidth from Kaplan and Sun (2017)
  # h.b <- ivqr.bw(p=tau, va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b.init=b2)
  h.b <- n^(-1/2)
  Lambda.derivative <- function(b){
    #State Variable
    g1 <- -Xt[,1:ncol(state)]+b[length(b)]*lX[,1:ncol(state)]
    #Free Variables
    g2 <- -Xt[,(ncol(state)+1):ncol(Xt)]+b[length(b)]*lX[,(ncol(state)+1):ncol(Xt)]
    # #Productivity
    g3 <- -(Lagphi-lX%*%b[1:(length(b)-1)])
    # # #Combined
    g <- cbind(g1, g2, g3, deparse.level = 0)
    return(g)
    }
  L.d <- Lambda.derivative(b2) 
  G.temp <- -t(array(data=Itilde.deriv(-Lambda(va=va, Xt=Xt, lX=lX, Lagphi=Lagphi, b=b2)/h.b),dim=dim(Z))*Z)%*%L.d
  G.b <- solve(G.temp)
  avar.b <- G.b%*%weight.mat%*%t(G.b)/n
  return(list(b2, as.numeric(diag(avar.b))))
}
#Load Chilean dataset
chile_panel <- read.csv('chile_panel.csv')
#Choose which industry to select
industries <- c(311, 381, 321, 331)
#Vector of quantiles
tau <- seq(0.1, 0.9, by=0.05)
#Store results for coefficients and standard errors
coefficients <- replicate(length(industries), list(array(0, dim=c(length(tau), 4))))
se <- replicate(length(industries), list(array(0, dim=c(length(tau), 4))))
#Estimate for different industries
#TO DO: Include Returns to Scale and test for returns to scale
for (ISIC in 1:length(industries)){
  print(industries[ISIC])
  chile <- subset(chile_panel, ciiu_3d==industries[ISIC])
  for (q in 1:length(tau)){
    print(tau[q])
    results <- QACF(tau=tau[q], va=chile$lnva, state=chile$lnk, free=cbind(chile$lnb, chile$lnw), proxy=chile$proxy_e, id=chile$id, time=chile$year, h=1e-6, b.init=c(0,0,0,0))
    coefficients[[ISIC]][q,] <- results[[1]]
    se[[ISIC]][q,] <- results[[2]]
    print(cbind(results[[1]], results[[2]]))
  }
}
#Prepare estimates for table in paper/presentation
require(xtable)
tau_table <- c(0.25, 0.5, 0.75)

estimates <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, coefficients), do.call(rbind, se))[,c(rbind(c(1:4), 4+(1:4)))]))
colnames(estimates) <- c('Tau','K',"se_K", 'Lw', "se_Lw", 'Lb', "se_Lb",'Rho', "se_Rho")

#Table Labels
ISIC_labels <- array(NA, length(tau_table)*length(industries)); ISIC_labels[seq(1, length(tau_table)*length(industries), by=length(tau_table))] <- industries
ISIC_labels[is.na(ISIC_labels)] <- ""

estimates_table <- subset(cbind(ISIC_labels, estimates[rep(tau, length(industries))%in%tau_table, ]), select=-c(Rho, se_Rho))
colnames(estimates_table) <- c("Industry (ISIC code)", "$\\tau$", "Coef.", 's.e.', "Coef.",'s.e.', "Coef.", 's.e.')

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4))
align(estimates_table) <- rep('c', 9)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Skilled Labor} & \\multicolumn{2}{c}{Unskilled Labor} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x)





