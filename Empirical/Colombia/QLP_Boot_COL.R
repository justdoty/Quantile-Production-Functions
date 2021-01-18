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
#Vector of quantiles of TFP
tfptau <- c(0.1, 0.25, 0.5, 0.9)
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tauvec[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#This gives the estimates using the unsampled data########################
#Elasticities from QLP
betahat <- array(0, dim=c(dZ, length(industries)))
#Elasticities from QR
qrhat <- array(0, dim=c(dZ, length(industries)))
#Difference between QLP and QR
qdifhat <- array(0, dim=c(dZ, length(industries)))
#Quantiles of TFP
QTFPhat <- array(0, dim=c(length(tfptau), length(industries)))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QLP
betaboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped Elasticities from QR
qrboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped differences between QLP and QR
qdifboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped quantiles of TFP
QTFPboot <- array(0, dim=c(R, length(tfptau), length(industries)))

for (isic in 1:length(industries)){
  COL <- filter(COLdata, str_detect(isic3, industries[isic]))
  soln <- QLP_Boot(tau=tau, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL, R=R, tfptau=tfptau)
  #Bootstrapped Estimates
  betaboot[,,isic] <- soln$betaboot
  qrboot[,,isic] <- soln$qrboot
  qdifboot[,,isic] <- soln$qdifboot
  QTFPboot[,,isic] <- soln$QTFPboot
  #Estimataes from Unsampled data
  betahat[,isic] <- soln$betahat
  qrhat[,isic] <- soln$qrhat
  qdifhat[,isic] <- soln$qdifhat
  QTFPhat[,isic] <- soln$QTFPhat
}
filename <- paste("PFQR/DATA/COL/QLP_Environments/QLP_Boot_COL_Q", id, ".RData", sep="")
save(tauvec, tfptau, dZ, betahat, qrhat, qdifhat, QTFPhat, betaboot, qrboot, qdifboot, QTFPboot, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) myjob.job
#Here length(tau)=19




