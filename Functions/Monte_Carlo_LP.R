# source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/Aux_Fun.R')
source('PFQR/FUN/Aux_Fun.R')
#For Paralelization
require(snow)
require(pracma)
#For First Stage QR
require(quantreg)
#Initialzie clusters
cl <- makeCluster(5)
################################################################################################
################################################################################################
################################################################################################
###############################DGP Parameters###################################################
#Specifications for Error Distributions
DGPs <- c("normal", "student", "log")
#MC Replications
nreps <- 1000
# nreps <- 3
#Vector of quantiles
tau <- seq(5,95,by=5)/100
# tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#Standard deviation of log wage process
siglnw <- 0
#Labor chosen at time timeb
timeb <- 0
#Standard deviation of optimization error
sigoptl <- 0.37
#Number of parameters to estimate in QGMM
dB <- 2
############################################################################################
###############Moment Equation Objective Function############################################
#############################################################################################
Lambda <- function(theta, mY, mX){
  resid <- mY-mX%*%theta[1:ncol(mX)]
  return(resid)
}
QLP <- function(theta, mY, mX, mZ, h, tau){
  xifit <- Lambda(theta=theta, mY=mY, mX=mX)
  gbar <- as.matrix(colMeans(mZ*repmat((Gfn(-xifit, h)-tau), 1, ncol(mZ))))
  W <- solve(tau*(1-tau)*t(mZ)%*%mZ/nrow(mZ))
  go <- nrow(mZ)*t(gbar)%*%W%*%gbar
  return(go)

}
############ LP Moment Equations############################
LP <- function(theta, mY, mX, mlX, mZ, fitphi, fitlagphi){
  A <- fitphi-mX%*%theta[1:ncol(mX)]
  B <- fitlagphi-mlX%*%theta[1:ncol(mX)]
  step1 <- lm(A~B)
  step1param <- as.numeric(coef(step1))
  wfit <- cbind(1,B)%*%step1param
  xsi <- mY-mX%*%theta[1:ncol(mX)]-wfit
  go <- crossprod(xsi)
  # mom <- mZ*array(xsi, dim(mZ))
  # momc <- colSums(mom^2)
  # go <- sum(momc)
  return(go)
}
##################################################################################
######################################################################################
#############################################################################################
#############################################################################################
####################Initialize Matrices to Store Results#####################################
# Number of Firms
n <- 1000
# Number of Time Periods
overallt <- 100
starttime <- 90
t <- overallt - starttime
#Store results for quantile estimators
resmat_LPQ <- array(0, dim=c(nreps, dB, length(tau), length(DGPs)))
resmat_QR <- array(0, dim=c(nreps, dB, length(tau), length(DGPs)))
resmat_TFP <- array(0, dim=c(nreps, n*t, length(DGPs)))
resmat_QTFP <- array(0, dim=c(nreps, n*t, length(DGPs)))
#Store results for LP estimator
resmat_LP <- array(0, dim=c(nreps, dB, length(DGPs)))
#Time entire code
overall.start.time <- Sys.time()
####################DGP########################################################
for (d in 1:length(DGPs)){
  for (j in 1:nreps){
    print(sprintf("DGP %i, Iteration %i", d, j))
    set.seed(123456+j)
    #Production Function Parameters
    alpha0 <- 0
    alphal <- 0.6
    alphak <- 0.4
    #Epsilons and omega ln(wage) process
    #For DGP 1
    sigeps <- 0.1
    #For DGP 2
    df <- 5
    siga <- sqrt((df-2)/df)*sigeps
    #For DGP3
    meanlog <- 0.15
    sdlog <- 0.1
    location <- log(meanlog^2/sqrt(meanlog^2+sdlog^2))
    shape <- sqrt(log(1+sdlog^2/meanlog^2))
    #Standard deviation of omega
    sigomg <- 0.3
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

    #Location Scale parameters
    eta0 <- 0
    etak <- c(0.6, 0.4, 0.5)
    etal <- c(-0.6, -0.4, -0.5)
  

    #Specification for Error Distribution for DGPs
    if (DGPs[d]=="normal"){
      epsdata <- matrix(rnorm(n*overallt, 0, sigeps), nrow=n, ncol=overallt)
      alphak0 <- alphak+etak[d]*qnorm(tau, 0, sigeps)
      alphal0 <- alphal+etal[d]*qnorm(tau, 0, sigeps)
      
    } else if (DGPs[d]=="student") {
      epsdata <- matrix(rt(n*overallt, df), nrow=n, ncol=overallt)
      epsdata <- siga*epsdata
      alphak0 <- alphak+etak[d]*qt(tau, 5)*siga
      alphal0 <- alphal+etal[d]*qt(tau, 5)*siga
    } else if (DGPs[d]=="log") {
      epsdata <- matrix(rlnorm(n*overallt, meanlog=location, sdlog=shape), nrow=n, ncol=overallt)
      alphak0 <- alphak+etak[d]*qlnorm(tau, meanlog=location, sdlog=shape)
      alphal0 <- alphal+etal[d]*qlnorm(tau, meanlog=location, sdlog=shape)
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
    #Specifies the form of heteroskedasticity
    het <- eta0+etal[d]*lnldata+etak[d]*lnkdata
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

    ##########################LP Estimation################################
    ################################################################################
    #First Stage########################################################
    firststage_LP <- lm(Output~Capital+Labor+Materials)
    phi0 <- firststage_LP$coefficients[1]
    LP_Labor <- firststage_LP$coefficients[3]
    resmat_LP[,,d][j,][2] <- LP_Labor
    phi <- fitted(firststage_LP)-as.matrix(Labor)%*%LP_Labor
    phi_LP <- phi
    dim(phi_LP) <- c(t, n)
    phi_Lag_1_LP <- c(phi_LP[1:(t-1),])
    phi_Con_LP <- c(phi_LP[2:t,])
    mY <- as.matrix(Output_Con-as.matrix(Labor_Con)%*%LP_Labor)
    #Matrix of Instruments
    mZ <- as.matrix(Capital_Con)
    #Contemporary Values
    mX <- as.matrix(Capital_Con)
    #Lag Values
    mlX <- as.matrix(Capital_Lag_1)
    results_LP <- optim(par=alphak, fn=function(theta) LP(theta, mY=mY, mX=mX, mlX=mlX, mZ=mZ, fitphi=phi_Con_LP, fitlagphi=phi_Lag_1_LP), gr=NULL, method="L-BFGS-B", lower=0, upper=1)$par
    #Estimated productivity
    wfit <- phi-Capital*results_LP
    resmat_TFP[,,d][j,] <- Output-cbind(Capital, Labor)%*%as.matrix(c(results_LP, LP_Labor))
    ############################################################
    resmat_LP[,,d][j,][1] <- results_LP[1]
    print("LP Estimates")
    print(resmat_LP[,,d][j,])
    ##################################Estimation############################################
    ##########################################################################################
    # ##########################################################################################
    clusterExport(cl, c('n','overallt','t','starttime','nreps', 'tau', 'dB', 'siglnw', 'timeb',
    'sigoptl', 'resmat_LPQ', 'rq', 'fitted',
    'Output', 'Capital', 'Labor', 'Materials',
    'Capital_Con','Capital_Lag_1', 'Labor_Lag_1', 'Labor_Con', 'Output_Con',
    'alphak', 'alphal', 'rho', 'j', 'd', 'DGPs', 'alphak0', 'alphal0', 'wfit', 'resmat_QR'), envir=environment())
    innerloop_QLP <- function(q){
      QR <- as.numeric(coef(rq(Output~Capital+Labor, tau=tau[q])))
      #Output
      mY <- Output-wfit
      #Matrix of Instruments
      mZ <- cbind(Capital, Labor)
      #Contemporary Values
      mX <- cbind(Capital, Labor)
      QRresults <- as.numeric(coef(rq(mY~mX-1, tau=tau[q])))
      QTFP <- Output-mX%*%QRresults
      #############################################################
      return(list(QRresults, QR[-1], QTFP))
      ##################################################################
    }
      #Optional for serial computing
      # resmat <- lapply(1:length(tau), innerloop_QLP)
      resmat <- parLapply(cl, 1:length(tau), innerloop_QLP)
      resmat_LPQ[,,,d][j,,] <- sapply(1:length(tau), function(l) rbind(resmat[[l]][[1]]))
      resmat_QR[,,,d][j,,] <- sapply(1:length(tau), function(l) rbind(resmat[[l]][[2]]))
      resmat_QTFP[,,d][j,] <- resmat[[(length(tau)+1)/2]][[3]]
      ####################################################################
      print("Q-GMM Estimates")
      print(t(resmat_LPQ[,,,d][j,,]))
  }
}
stopCluster(cl); print("Cluster stopped.")
print(Sys.time()-overall.start.time)
# #Store True Values in Rdata environment
#DGP 1
betak1 <- alphak+etak[1]*qnorm(tau, 0, sigeps); betal1 <- alphal+etal[1]*qnorm(tau, 0, sigeps)
beta1 <- cbind(betak1, betal1)
#DGP 2
betak2 <- alphak+etak[2]*qt(tau, 5)*siga; betal2 <- alphal+etal[2]*qt(tau, 5)*siga
beta2 <- cbind(betak2, betal2)
#DGP 3
betak3 <- alphak+etak[3]*qlnorm(tau, meanlog=location, sdlog=shape); betal3 <- alphal+etal[3]*qlnorm(tau, meanlog=location, sdlog=shape)
beta3 <- cbind(betak3, betal3)
#Combined
beta <- rbind(beta1, beta2, beta3)
#Calculate Estimates for Mean TFP densities and Confidence Intervals
#QTFP
QTFP <- array(0, dim=c(n*t, length(DGPs)))
QTFPUB <- array(0, dim=c(n*t, length(DGPs)))
QTFPLB <- array(0, dim=c(n*t, length(DGPs)))
#TFP
TFP <- array(0, dim=c(n*t, length(DGPs)))
TFPUB <- array(0, dim=c(n*t, length(DGPs)))
TFPLB <- array(0, dim=c(n*t, length(DGPs)))
for (i in 1:length(DGPs)){
  #QTFP
  QTFP[,i] <- colMeans(resmat_QTFP[,,i])
  QTFPUB[,i] <- apply(resmat_QTFP[,,i], 2, function(q) quantile(q, probs=0.975))
  QTFPLB[,i] <- apply(resmat_QTFP[,,i], 2, function(q) quantile(q, probs=0.025))
  #TFP
  TFP[,i] <- colMeans(resmat_TFP[,,i])
  TFPUB[,i] <- apply(resmat_TFP[,,i], 2, function(q) quantile(q, probs=0.975))
  TFPLB[,i] <- apply(resmat_TFP[,,i], 2, function(q) quantile(q, probs=0.025))
}
#Save Results
save(nreps, DGPs, resmat_LP, resmat_LPQ, resmat_QR, QTFP, QTFPUB, QTFPLB, TFP, TFPUB, TFPLB, beta, tau, dB, file="PFQR/SIM/simulation_LP.Rdata")
















