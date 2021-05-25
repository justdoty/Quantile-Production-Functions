#This the file that calls the functions to estimate the production function for the COL data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QACF_Boot.R')
require(stringr)
#Load COL dataset
COLdata <- read.csv("PFQR/DATA/COL/GNRCOLdata.csv") %>% transmute(id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), isic3=isic3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
#Vector of quantiles
tauvec <- seq(5, 95, length.out=19)/100
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
isic <- industries[id]
#The number of bootstrap replications to be used in QACF defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#This gives the estimates using the unsampled data########################
#Elasticities from QACF
QACFbetahat <- array(0, dim=c(dZ, tau_n))
ACFhat <- array(0, dim=c(dZ, tau_n))
#Elasticities from QR
QACFqrhat <- array(0, dim=c(dZ, tau_n))
#Difference between QACF and QR
QACFqdifhat <- array(0, dim=c(dZ, tau_n))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QACF
QACFbetaboot <- array(0, dim=c(R, dZ, tau_n))
ACFboot <- array(0, dim=c(R, dZ, tau_n))
#Bootstrapped Elasticities from QR
QACFqrboot <- array(0, dim=c(R, dZ, tau_n))
#Bootstrapped differences between QACF and QR
QACFqdifboot <- array(0, dim=c(R, dZ, tau_n))

for (j in 1:tau_n){
  print(sprintf("tau=%s", tauvec[j]))
  COL <- filter(COLdata, str_detect(isic3, isic))
  soln <- QACF_Boot(tau=tauvec[j], idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL, R=R)
  #Bootstrapped Estimates
  QACFbetaboot[,,j] <- soln$betaboot
  ACFboot[,,j] <- soln$ACFboot
  QACFqrboot[,,j] <- soln$qrboot
  QACFqdifboot[,,j] <- soln$qdifboot
  #Estimataes from Unsampled data
  QACFbetahat[,j] <- soln$betahat
  ACFhat[,j] <- soln$ACFhat
  QACFqrhat[,j] <- soln$qrhat
  QACFqdifhat[,j] <- soln$qdifhat
}
filename <- paste("PFQR/DATA/COL/QACF_Environments/QACF_Boot_COL_ISIC", id, ".RData", sep="")
save(tauvec, dZ, QACFbetahat, ACFhat, QACFqrhat, QACFqdifhat, QACFbetaboot, ACFboot, QACFqrboot, QACFqdifboot, file=filename)
###############################################################################################
#Repeat for QLP for Gross-Output Production Function
#################################################################################################
source('PFQR/FUN/QLP_Boot.R')
dZLP <- 3
#This gives the estimates using the unsampled data########################
#Elasticities from QLP
QLPbetahat <- array(0, dim=c(dZLP, tau_n))
LPhat <- array(0, dim=c(dZLP, tau_n))
#Elasticities from QR
QLPqrhat <- array(0, dim=c(dZLP, tau_n))
#Difference between QLP and QR
QLPqdifhat <- array(0, dim=c(dZLP, tau_n))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QLP
QLPbetaboot <- array(0, dim=c(R, dZLP, tau_n))
LPboot <- array(0, dim=c(R, dZLP, tau_n))
#Bootstrapped Elasticities from QR
QLPqrboot <- array(0, dim=c(R, dZLP, tau_n))
#Bootstrapped differences between QLP and QR
QLPqdifboot <- array(0, dim=c(R, dZLP, tau_n))

for (j in 1:tau_n){
  print(sprintf("tau=%s", tauvec[j]))
  COL <- filter(COLdata, str_detect(isic3, isic))
  soln <- QLP_Boot(tau=tauvec[j], idvar=COL$id, timevar=COL$year, Y=COL$lny, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZLP, binit=NULL, R=R)
  #Bootstrapped Estimates
  QLPbetaboot[,,j] <- soln$betaboot
  LPboot[,,j] <- soln$LPboot
  QLPqrboot[,,j] <- soln$qrboot
  QLPqdifboot[,,j] <- soln$qdifboot
  #Estimataes from Unsampled data
  QLPbetahat[,j] <- soln$betahat
  LPhat[,j] <- soln$LPhat
  QLPqrhat[,j] <- soln$qrhat
  QLPqdifhat[,j] <- soln$qdifhat
}
filename <- paste("PFQR/DATA/COL/QLP_Environments/QLP_Boot_COL_ISIC", id, ".RData", sep="")
save(tauvec, dZLP, QLPbetahat, LPhat, QLPqrhat, QLPqdifhat, QLPbetaboot, LPboot, QLPqrboot, QLPqdifboot, file=filename)
#HPC Job Submissions for batches: qsub -t 1:length(industries) QACF_Boot_COL.job
#Here length(industries)=4


