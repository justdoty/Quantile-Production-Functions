#This the file that calls the functions to estimate the production function for the COL data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QACF_Boot.R')
require(stringr)
#Load COL dataset
COLdata <- read.csv("PFQR/DATA/COL/COLdata.csv")
COLdata <- COLdata %>% group_by(id) %>% filter(n()>=3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
#Vector of quantiles
tauvec <- seq(5, 95, length.out=19)/100
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tauvec[id]
#The number of bootstrap replications to be used in QACF defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#This gives the estimates using the unsampled data########################
#Elasticities from QACF
QACFbetahat <- array(0, dim=c(dZ, length(industries)))
ACFhat <- array(0, dim=c(dZ, length(industries)))
#Elasticities from QR
QACFqrhat <- array(0, dim=c(dZ, length(industries)))
#Difference between QACF and QR
QACFqdifhat <- array(0, dim=c(dZ, length(industries)))
#QTFP Estimates
QACFTFPhat <- list()
#TFP Estimates
ACFTFPhat <- list()
#Omega Estimates
ACFomegahat <- list()
#Ex-post Shock Estiamtes
ACFexpost <- list()
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QACF
QACFbetaboot <- array(0, dim=c(R, dZ, length(industries)))
ACFboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped Elasticities from QR
QACFqrboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped differences between QACF and QR
QACFqdifboot <- array(0, dim=c(R, dZ, length(industries)))

for (isic in 1:length(industries)){
  COL <- filter(COLdata, str_detect(isic3, industries[isic]))
  soln <- QACF_Boot(tau=tau, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL, R=R)
  #Bootstrapped Estimates
  QACFbetaboot[,,isic] <- soln$betaboot
  ACFboot[,,isic] <- soln$ACFboot
  QACFqrboot[,,isic] <- soln$qrboot
  QACFqdifboot[,,isic] <- soln$qdifboot
  #Estimataes from Unsampled data
  QACFbetahat[,isic] <- soln$betahat
  ACFhat[,isic] <- soln$ACFhat
  QACFqrhat[,isic] <- soln$qrhat
  QACFqdifhat[,isic] <- soln$qdifhat
  QACFTFPhat[[isic]] <- soln$QTFPhat
  ACFTFPhat[[isic]] <- soln$TFPhat
  ACFomegahat[[isic]] <- soln$omegahat
  ACFexpost[[isic]] <- soln$expost
}
filename <- paste("PFQR/DATA/COL/QACF/QACF_Environments/QACF_Boot_COL_Q", id, ".RData", sep="")
save(tauvec, dZ, QACFbetahat, ACFhat, QACFqrhat, QACFqdifhat, QACFTFPhat, ACFTFPhat, ACFomegahat, QACFbetaboot, ACFboot, QACFqrboot, QACFqdifboot, ACFexpost, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) myjob.job
#Here length(tau)=19



