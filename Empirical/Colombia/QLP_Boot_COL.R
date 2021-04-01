#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP_Boot.R')
require(stringr)
#Load US dataset
COLdata <- read.csv("PFQR/DATA/COL/COLdata.csv")
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
#Vector of quantiles
tauvec <- seq(5, 95, length.out=19)/100
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tauvec[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#This gives the estimates using the unsampled data########################
#Elasticities from QLP
QLPbetahat <- array(0, dim=c(dZ, length(industries)))
LPhat <- array(0, dim=c(dZ, length(industries)))
#Elasticities from QR
QLPqrhat <- array(0, dim=c(dZ, length(industries)))
#Difference between QLP and QR
QLPqdifhat <- array(0, dim=c(dZ, length(industries)))
#QTFP Estimates
QLPTFPhat <- list()
#TFP Estimates
LPTFPhat <- list()
#Omega Estimates
LPomegahat <- list()
#Expost shock estimates
LPexpost <- list()
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QLP
QLPbetaboot <- array(0, dim=c(R, dZ, length(industries)))
LPboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped Elasticities from QR
QLPqrboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped differences between QLP and QR
QLPqdifboot <- array(0, dim=c(R, dZ, length(industries)))

for (isic in 1:length(industries)){
  COL <- filter(COLdata, str_detect(isic3, industries[isic]))
  soln <- QLP_Boot(tau=tau, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL, R=R)
  #Bootstrapped Estimates
  QLPbetaboot[,,isic] <- soln$betaboot
  LPboot[,,isic] <- soln$LPboot
  QLPqrboot[,,isic] <- soln$qrboot
  QLPqdifboot[,,isic] <- soln$qdifboot
  #Estimataes from Unsampled data
  QLPbetahat[,isic] <- soln$betahat
  LPhat[,isic] <- soln$LPhat
  QLPqrhat[,isic] <- soln$qrhat
  QLPqdifhat[,isic] <- soln$qdifhat
  QLPTFPhat[[isic]] <- soln$QTFPhat
  LPTFPhat[[isic]] <- soln$TFPhat
  LPomegahat[[isic]] <- soln$omegahat
  LPexpost[[isic]] <- soln$expost
}
filename <- paste("PFQR/DATA/COL/QLP_Environments/QLP_Boot_COL_Q", id, ".RData", sep="")
save(tauvec, dZ, QLPbetahat, LPhat, QLPqrhat, QLPqdifhat, QLPTFPhat, LPTFPhat, LPomegahat, QLPbetaboot, LPboot, QLPqrboot, QLPqdifboot, LPexpost, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) myjob.job
#Here length(tau)=19




