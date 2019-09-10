setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
source('gmmq.R')
source('ivqr_gmm.R')
#For Paralelization
require(snow)
#For MM 
require(pracma)
#For First Stage QR
require(quantreg)
#For GMM optimiztion
require(GenSA)
#Initialzie clusters
cl <- makeCluster(4)
################################################################################################
################################################################################################
################################################################################################
###############################DGP Parameters###################################################
#Specifications for Error Distributions
DGPs <- c("normal", "log-normal")
#MC Replications
nreps <- 2
#Vector of quantiles
tau <- seq(0.25, 0.75, by=0.25)
#Standard deviation of log wage process
siglnw <- 0.1
#Labor chosen at time timeb
timeb <- 0.5
#Standard deviation of optimization error
sigoptl <- 0
#Number of parameters to estimate
dB <- 3
###############Moment Equation for Q-GMM#######################################################
Moment <- function(x,b){
      Moment <- x[,1]-b[1]*x[,2]-b[2]*x[,3]-b[3]*(x[,4]-b[1]*x[,5]-b[2]*x[,6])
      return(Moment)
    }
#############Jacobian of Moment Equation for Q-GMM##################################
Moment.der <- function(x,b){
      #Capital
      g1 <- -x[,2]+b[3]*x[,5]
      #Labor
      g2 <- -x[,3]+b[3]*x[,6]
      #Productivity
      g3 <- -(x[,4]-b[1]*x[,5]-b[2]*x[,6])
      #Combined
      g <- cbind(g1, g2, g3, deparse.level = 0)
      return(g)
    }
#############ACF Objective Function for Exact Identification########################
ACF_MM <- function(x,z,b){
  omega2 <- x[,1]-x[,2]*b[1]-x[,3]*b[2]
  omega1 <- x[,4]-x[,5]*b[1]-x[,6]*b[2]
  omega.fit <- lm(omega2~omega1)
  Moment <- as.matrix(omega2-fitted(omega.fit))
  Obj <- colMeans(z*array(data=-Moment, dim=dim(z)))
  return(Obj)
} 
#############ACF Moment Equations for Over-Identification############################
ACF_GMM <- function(x,z,b){
  omega2 <- x[,1]-x[,2]*b[1]-x[,3]*b[2]
  omega1 <- x[,4]-x[,5]*b[1]-x[,6]*b[2]
  omega.fit <- lm(omega2~omega1-1)
  Moment <- omega2-fitted(omega.fit)
  Obj <- z*array(data=-Moment, dim=dim(z))
  return(Obj)
}  
##################################################################################
######################################################################################
#############################################################################################
#############################################################################################
####################Initialize Matrices to Store Results#####################################
resmat <- array(0, dim=c(nreps, dB, length(tau), length(DGPs)))
resmat_ACF <- array(0, dim=c(nreps, 2, length(DGPs)))
#Time entire code
overall.start.time <- Sys.time()
####################DGP########################################################
for (d in 1:length(DGPs)){
  for (j in 1:nreps){
    print(sprintf("DGP %i, Iteration %i", d, j))
    set.seed(j)
    # Number of Firms
    n <- 1000
    # Number of Time Periods
    overallt <- 100
    starttime <- 90
    t <- overallt - starttime
    #Production Function Parameters
    alpha0 <- 0
    alphal <- 0.6
    alphak <- 0.4
    #Epsilons and omega ln(wage) process
    sigeps <- 0.1
    mlog <- 0.4; slog <- 0.2
    meanlog <- log(mlog^2/sqrt(mlog^2+slog^2))
    sdlog <- sqrt(log(1+(slog^2)/mlog^2))
    siglog <- 0.2
    sigomg <- 0.3 #standard deviation of omega
    rho <- 0.7 #AR(1) coefficient for omega
    sigxi <- sqrt((1-rho^2)*sigomg^2)
    rholnw <- 0.3 # AR(1) coefficient for ln(wage)
    #Matrices to store data
    lnkdata <- matrix(0, n, overallt) #ln(capital)
    lnldata <- matrix(0, n, overallt) #ln(labor)
    lnmdata <- matrix(0, n, overallt) #ln(intermediate input)
    lnwdata <- matrix(0, n, overallt) #ln(wage)
    lnpdata <- matrix(0, n, overallt) #ln(output price)
    lnydata <- matrix(0, n, overallt) #ln(output)
    omgdata <- matrix(0, n, overallt) #omega(t)
    omgdataminusb <- matrix(0, n, overallt) #omega(t-b)

    #Location Scale Parameters
    etak <- 0.7
    etal <- 0.6
    etaomega <- 0.1

    #Specification for Error Distribution for DGPs
    if (DGPs[d]=="normal"){
      epsdata <- matrix(rnorm(n*overallt, 0, sigeps), nrow=n, ncol=overallt)
      sigxidata <- matrix(rnorm(n*overallt, 0, sigxi), nrow=n, ncol=overallt)
      alphak0 <- alphak+etak*qnorm(tau, 0, sigeps)
      alphal0 <- alphal+etal*qnorm(tau, 0, sigeps)
      
    } else if (DGPs[d]=="log-normal") {
      epsdata <- matrix(rlnorm(n*overallt, meanlog, sdlog), nrow=n, ncol=overallt)
      sigxidata <- matrix(rnorm(n*overallt, 0, sigxi), nrow=n, ncol=overallt)
      alphak0 <- alphak+etak*qlnorm(tau, meanlog, sdlog)
      alphal0 <- alphal+etal*qlnorm(tau, meanlog, sdlog)
    } else {
      print("Error: Unspecified DGP")
    }


    #subdividing the AR(1) process
    rhofirst <- rho^(1-timeb)
    rhosecond <- rho^(timeb)
    sigxifirst <- sqrt((1-rhofirst^2)*sigomg^2)
    sigxisecond <- sqrt((1-rhosecond^2)*sigomg^2) #Standard deviation of innovation in omega

    sigxilnw <- sqrt((1-rholnw^2)*siglnw^2) #standard deviation of innovation in lnw(wage)

    #Period 0 values of omega and ln(wage)
    omgdata0 <- matrix(rnorm(n,0,sigomg),nrow=n,ncol=1)
    lnwdata0 <- matrix(rnorm(n,0,siglnw),nrow=n,ncol=1)

    #Period 1-b values of omega and period 1 values of omega and ln(wage)
    omgdataminusb[,1] <- rhofirst*omgdata0+matrix(rnorm(n,0,sigxifirst),nrow=n,ncol=1)
    omgdata[,1] <- rhosecond*omgdataminusb[,1]+matrix(rnorm(n,0,sigxisecond),nrow=n,ncol=1)
    lnwdata[,1] <- rholnw*lnwdata0 + matrix(rnorm(n,0,sigxilnw),nrow=n,ncol=1)


    #Simulate values of omega and ln(wage) for rest of time periods
    for (s in 2:overallt){
      omgdataminusb[,s] <- rhofirst*omgdata[,s-1] + matrix(rnorm(n,0,sigxifirst),nrow=n,ncol=1)
      omgdata[,s] <- rhosecond*omgdataminusb[,s] + matrix(rnorm(n,0,sigxisecond),nrow=n,ncol=1)
      lnwdata[,s] <- rholnw*lnwdata[,s-1] + matrix(rnorm(n,0,sigxilnw),nrow=n,ncol=1)
    }


    #Intital ln(capital) level (close to 0)
    lnkdata[,1] <- matrix(-100,n,1)
    #Discount Rate for DP problem
    disc <- 0.95
    #Depreciation rate of capital
    delta <- 0.2
    #Variation in capital adjustment costs across firms
    sigb <- 0.6
    #See page 36 ACF 1/Phi(i) is distributed lognormally 
    #across firms but constant overtime with sd 0.6
    oneoverbiadj <- exp(rnorm(n,0,sigb))

    #Simplifying components of the optimal investment rule
    #Square bracket component
    squarebracketterm <- (alphal^(alphal/(1-alphal)))*exp(0.5*alphal^2*sigoptl^2) -
    (alphal^(1/(1-alphal)))*exp(0.5*sigoptl^2)
    #Constant term in front of sum including squarebracketterm
    const1 <- disc*(alphak/(1-alphal))*(exp(alpha0)^(1/(1-alphal)))*
    squarebracketterm
    vec1 <- (disc*(1-delta))^seq(100)
    vec2 <- cumsum(rholnw^(2*seq(100)))
    vec3 <- (sigxi^2)*cumsum(rho^(2*seq(0,99,1)))
    vec3 <- cbind(sigxi^2 * 0,cumsum(rho^(2*(seq(100)-1))) )
    expterm3 <- exp(0.5*((-alphal)/(1-alphal))^2*((sigxilnw^2)*vec2))
    expterm4 <- exp(0.5*(1/(1-alphal))^2*rhosecond^2*
    ((sigxifirst^2)*rho^(2*seq(100))+vec3))
    expterm5 <- exp((1/(1-alphal))*(1/2)*sigxisecond^2)
    #Compute Optimal Investment and Capital stock for all firms over time
    investmat <- matrix(NA, n, overallt)

    for (i in 1:n){
        for (s in 1:overallt){
          expterm1 <- exp((1/(1-alphal))*omgdata[i,s]*rho^(seq(100)))
          expterm2 <- exp(((-alphal)/(1-alphal))*lnwdata[i,s]*(rholnw^seq(100)))
          investmat[i,s] <- oneoverbiadj[i]*const1*expterm5*sum(vec1*expterm1*
            expterm2*expterm3*expterm4)

          if (s >= 2){
            lnkdata[i,s] <- log((1-delta)*exp(lnkdata[i,s-1])+
              (1-0*runif(1))*investmat[i,s-1])
          }
        }
      }


    #Generate levels of labor input
    for (s in 1:overallt){
      lnldata[,s] <- ((sigxisecond^2)/2+log(alphal)+alpha0+
        rhosecond*omgdataminusb[,s]-lnwdata[,s]+lnpdata[,s]+
        (alphak)*lnkdata[,s])/(1-alphal)
    }

    #Potential Optimization Error
    truelnldata <- lnldata
    lnldata <- lnldata + matrix(rnorm(n*overallt,0,sigoptl),n,overallt)

    #Output and Materials
    het <- etal*lnldata+etak*lnkdata+etaomega*omgdata
    lnydata <- alpha0 + alphal*lnldata + alphak*lnkdata + omgdata + het*epsdata
    lnmdata <- alpha0 + alphal*truelnldata + alphak*lnkdata + omgdata

    #Stack data across firms (all the data)
    Capital <- c(t(lnkdata[,(starttime+1):overallt]))
    Labor <- c(t(lnldata[,(starttime+1):overallt]))
    Materials <- c(t(lnmdata[,(starttime+1):overallt]))
    Wage <- c(t(lnwdata[,(starttime+1):overallt]))
    Price <- c(t(lnpdata[,(starttime+1):overallt]))
    Output <- c(t(lnydata[,(starttime+1):overallt]))
    Productivity <- c(t(omgdata[,(starttime+1):overallt]))
    Productivity_t_minus_b <- c(t(omgdataminusb[,(starttime+1):overallt]))
    True_Labor <- c(t(truelnldata[,(starttime+1):overallt]))
    Epsilon <- c(t(epsdata[,(starttime+1):overallt]))

    #Stack data across firms (lagged data)
    Capital_Lag_1 <- c(t(lnkdata[,(starttime+1):(overallt-1)]))
    Labor_Lag_1 <- c(t(lnldata[,(starttime+1):(overallt-1)]))
    Labor_Lag_2 <- c(t(lnldata[,(starttime):(overallt-2)]))
    Materials_Lag_1 <- c(t(lnmdata[,(starttime+1):(overallt-1)]))
    Wage_Lag_1 <- c(t(lnwdata[,(starttime+1):(overallt-1)]))
    Price_Lag_1 <- c(t(lnpdata[,(starttime+1):(overallt-1)]))
    Output_Lag_1 <- c(t(lnydata[,(starttime+1):(overallt-1)]))
    Productivity_Lag_1 <- c(t(omgdata[,(starttime+1):(overallt-1)]))
    Productivity_t_minus_b_Lag_1 <- c(t(omgdataminusb[,(starttime+1):(overallt-1)]))
    True_Labor_Lag_1 <- c(t(truelnldata[,(starttime+1):(overallt-1)]))

    #Stack data across firms (contemporaneous data)
    Capital_Con <- c(t(lnkdata[,(starttime+2):(overallt)]))
    Labor_Con <- c(t(lnldata[,(starttime+2):(overallt)]))
    Materials_Con <- c(t(lnmdata[,(starttime+2):(overallt)]))
    Wage_Con <- c(t(lnwdata[,(starttime+2):(overallt)]))
    Price_Con <- c(t(lnpdata[,(starttime+2):(overallt)]))
    Output_Con <- c(t(lnydata[,(starttime+2):(overallt)]))
    Productivity_Con <- c(t(omgdata[,(starttime+2):(overallt)]))
    Productivity_t_minus_b_Con <- c(t(omgdataminusb[,(starttime+2):(overallt)]))
    True_Labor_Con <- c(t(truelnldata[,(starttime+2):(overallt)]))

    ##################################ACF Estimation################################
    ################################################################################
    ###############First Stage########################################################
    firststage_ACF <- lm(Output~Capital+Labor+Materials)
    phiacf_ACF <- fitted(firststage_ACF)
    dim(phiacf_ACF) <- c(t, n)
    phiacf_Lag_1_ACF <- c(phiacf_ACF[1:(t-1),])
    phiacf_Con_ACF <- c(phiacf_ACF[2:t,])
    ACF_Z <- cbind(1, Capital_Con, Capital_Lag_1, Labor_Lag_1, Labor_Lag_2)
    ACF_X <- cbind(phiacf_Con_ACF, Capital_Con, Labor_Con, phiacf_Lag_1_ACF, Capital_Lag_1, Labor_Lag_1)
    #Method of Moments - Exact Identification
    # results_ACF <- newtonsys(Ffun=function(b) ACF_MM(x=ACF_X, z=ACF_Z, b), 
    #   Jfun=NULL, x0=c(alphak, alphal), maxiter=400)$zero
    # Quantile GMM - Exact Identification and Over-Identification
    obj.fn <- function(b){
      momi <- ACF_GMM(x=ACF_X, z=ACF_Z, b)
      return(nrow(momi)*colMeans(momi)%*%inv(var(momi))%*%as.matrix(colMeans(momi)))
    }
    results_ACF <- GenSA(par=c(alphak, alphal), fn=obj.fn, lower=c(0,0),
      upper=c(1, 1), control=list(max.time=5))$par
    ############################################################
    resmat_ACF[,,d][j,] <- results_ACF
    print("ACF Estimates")
    print(resmat_ACF[,,d][j,])

    ##################################Estimation############################################
    ##########################################################################################
    # ##########################################################################################
    clusterExport(cl, c('n','overallt','t','starttime','nreps', 'tau', 'dB', 'siglnw', 'timeb',
    'sigoptl', 'gmmq', 'Itilde.KS17',
    'Itilde.deriv.KS17', 'resmat', 'Moment', 'Moment.der', 'rq', 'fitted',
    'Output', 'Capital', 'Labor', 'Materials',
    'Capital_Con','Capital_Lag_1', 'Labor_Lag_1', 'Labor_Lag_2', 'Labor_Con', 'Output_Con',
    'alphak', 'alphal', 'rho', 'j', 'd', 'DGPs', 'alphak0', 'alphal0',
    'ivqr.gmm', 'LRV.est.fn', 'uniform.fn', 'IDENTITY.FLAG', 'ZZ.FLAG',
    'GenSA', 'repmat'), envir=environment())
    innerloop <- function(q){
      firststage <- rq(Output~Capital+Labor+Materials, tau=tau[q])
      phiacf <- fitted(firststage)
      dim(phiacf) <- c(t, n)
      phiacf_Lag_1 <- c(phiacf[1:(t-1),])
      phiacf_Con <- c(phiacf[2:t,])
      Z <- cbind(Capital_Con, Labor_Lag_1, phiacf_Lag_1)
      X <- cbind(Output_Con, Capital_Con, Labor_Con, phiacf_Lag_1, Capital_Lag_1, Labor_Lag_1, phiacf_Con)
      # #Method of Moments - Exact Identification
      # results <- gmmq(tau=tau[q], X=X, Z=Z, dB=dB, Lambda=Moment,
      #   Lambda.derivative=Moment.der, h=0, VERBOSE=FALSE, RETURN.Z=FALSE,
      #   b.init=c(alphak0[q], alphal0[q], rho))
      #Quantile GMM - Exact Identification and Over-Identification
      results <- ivqr.gmm(tau=tau[q], X=X, Z=Z, h=0.001, dB=dB, max.time=5, upper=c(1, 1, 1),
        lower=c(0, 0, 0), structure='iid',
        LRV.kernel='uniform', Lambda=Moment, Lambda.derivative=Moment.der, b.init=c(alphak0[q], alphal0[q], rho))
      #############################################################
      resmat[,,,d][,,q][j,] <- t(results$b)
      return(resmat[,,,d][,,q][j,])
      ##################################################################
    }
    #Optional for serial computing
    # resmat[,,,d][j,,] <- matrix(unlist(lapply(1:length(tau), innerloop)), nrow=length(tau), ncol=dB)
    q.time <- proc.time()
    resmat[,,,d][j,,] <- matrix(unlist(parLapply(cl, 1:length(tau), innerloop)), nrow=length(tau), ncol=dB)
    ####################################################################
    print("Q-GMM Estimates")
    print(t(resmat[,,,d][j,,]))
    print(proc.time()-q.time)
  }
}
stopCluster(cl); print("Cluster stopped.")
print(Sys.time()-overall.start.time)
#Store True Values in Rdata environment
alphak1 <- alphak+etak*qnorm(tau, 0, sigeps); alphal1 <- alphal+etal*qnorm(tau, 0, sigeps)
alpha1 <- cbind(alphak1, alphal1)
alphak2 <- alphak+etak*qlnorm(tau, meanlog, sdlog); alphal2 <- alphal+etal*qlnorm(tau, meanlog, sdlog)
alpha2 <- cbind(alphak2, alphal2)
alpha <- rbind(alpha1, alpha2)
#Save Results
# save(nreps, DGPs, resmat, resmat_ACF, alpha, tau, dB, file="simulation_data.Rdata")

























