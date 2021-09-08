#This the file that calls the functions to estimate the production function for the COL data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QACF_Boot.R')
require(stringr)
#Load COL dataset
COL_panel <- read.csv("PFQR/DATA/COL/COLdata.csv")
COLdata <- transmute(COL_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), ExB=ExB, ImB=ImB, AdvB=AdvB, ExC=ExC, ImC=ImC, AdvC=AdvC, isic3=isic3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "381", "321", All)
#Vector of quantiles
tauvec <- seq(5, 95, length.out=19)/100
tau_n <- length(tauvec)
id <- as.numeric(commandArgs(TRUE)[1])
isic <- industries[id]
#Filter COL data per industry
COL <- filter(COLdata, str_detect(isic3, isic))
#The number of bootstrap replications to be used in QACF defined below
R <- 500
# The number of parameters being estimated
dZ <- 2
# This gives the estimates using the unsampled data########################
# Elasticities from QACF
QACFbetahat <- array(0, dim=c(dZ, tau_n))
ACFhat <- array(0, dim=c(dZ, tau_n))
#Elasticities from QR
QACFqrhat <- array(0, dim=c(dZ, tau_n))
#Difference between QACF and QR
QACFqdifhat <- array(0, dim=c(dZ, tau_n))
#Productivity Differentials
QACFDSPB <- array(0, dim=c(3, tau_n))
QACFDSPC <- array(0, dim=c(3, tau_n))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QACF
QACFbetaboot <- array(0, dim=c(R, dZ, tau_n))
ACFboot <- array(0, dim=c(R, dZ, tau_n))
#Bootstrapped Elasticities from QR
QACFqrboot <- array(0, dim=c(R, dZ, tau_n))
#Bootstrapped differences between QACF and QR
QACFqdifboot <- array(0, dim=c(R, dZ, tau_n))
#Bootstrapped Productivity Differentials
QACFDSPB_boot <- array(0, dim=c(R, 3, tau_n))
QACFDSPC_boot <- array(0, dim=c(R, 3, tau_n))
for (j in 1:tau_n){
  print(sprintf("tau=%s", tauvec[j]))
  soln <- QACF_Boot(tau=tauvec[j], idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL, R=R, XC=cbind(COL$ExC, COL$ImC, COL$AdvC), XB=cbind(COL$ExB, COL$ImB, COL$AdvB))
  #Bootstrapped Estimates
  QACFbetaboot[,,j] <- soln$betaboot
  QACFqrboot[,,j] <- soln$qrboot
  QACFqdifboot[,,j] <- soln$qdifboot
  #Estimates from Productivity Differentials
  QACFDSPB_boot[,,j] <- soln$DSPBboot
  QACFDSPC_boot[,,j] <- soln$DSPCboot
  #Estimataes from Unsampled data
  QACFbetahat[,j] <- soln$betahat
  QACFqrhat[,j] <- soln$qrhat
  QACFqdifhat[,j] <- soln$qdifhat
  #Estimates from Productivity Differentials
  QACFDSPB[,j] <- soln$DSPB
  QACFDSPC[,j] <- soln$DSPC
  if (j==1){
    ACFhat <- soln$ACFhat
    ACFboot <- soln$ACFboot
    #Estimates from Productivity Differentials
    ACFPB <- soln$ACFPB
    ACFPC <- soln$ACFPC
    ACFPB_boot <- soln$ACFPBboot
    ACFPC_boot <- soln$ACFPCboot
  }
}
filename <- paste("PFQR/DATA/COL/QACF_Environments/QACF_Boot_COL_ISIC", id, ".RData", sep="")
save(tauvec, dZ, QACFbetahat, ACFhat, QACFqrhat, QACFqdifhat, QACFbetaboot, 
  ACFboot, QACFqrboot, QACFqdifboot, QACFDSPB, QACFDSPC, 
  QACFDSPB_boot, QACFDSPC_boot, ACFPB, ACFPC, ACFPB_boot, ACFPC_boot, file=filename)
###############################################################################################
#Repeat for QLP for Gross-Output Production Function
#################################################################################################
source('PFQR/FUN/QLP_Boot.R')
dZLP <- 3
#This gives the estimates using the unsampled data########################
#Elasticities from QACF
QLPbetahat <- array(0, dim=c(dZLP, tau_n))
LPhat <- array(0, dim=c(dZLP, tau_n))
#Elasticities from QR
QLPqrhat <- array(0, dim=c(dZLP, tau_n))
#Difference between QLP and QR
QLPqdifhat <- array(0, dim=c(dZLP, tau_n))
#Productivity Differentials
QLPDSPB <- array(0, dim=c(3, tau_n))
QLPDSPC <- array(0, dim=c(3, tau_n))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QLP
QLPbetaboot <- array(0, dim=c(R, dZLP, tau_n))
LPboot <- array(0, dim=c(R, dZLP, tau_n))
#Bootstrapped Elasticities from QR
QLPqrboot <- array(0, dim=c(R, dZLP, tau_n))
#Bootstrapped differences between QLP and QR
QLPqdifboot <- array(0, dim=c(R, dZLP, tau_n))
#Bootstrapped Productivity Differentials
QLPDSPB_boot <- array(0, dim=c(R, 3, tau_n))
QLPDSPC_boot <- array(0, dim=c(R, 3, tau_n))
for (j in 1:tau_n){
  print(sprintf("tau=%s", tauvec[j]))
  soln <- QLP_Boot(tau=tauvec[j], idvar=COL$id, timevar=COL$year, Y=COL$lny, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZLP, binit=NULL, R=R, XC=cbind(COL$ExC, COL$ImC, COL$AdvC), XB=cbind(COL$ExB, COL$ImB, COL$AdvB))
  #Bootstrapped Estimates
  QLPbetaboot[,,j] <- soln$betaboot
  QLPqrboot[,,j] <- soln$qrboot
  QLPqdifboot[,,j] <- soln$qdifboot
  #Estimates from Productivity Differentials
  QLPDSPB_boot[,,j] <- soln$DSPBboot
  QLPDSPC_boot[,,j] <- soln$DSPCboot
  #Estimataes from Unsampled data
  QLPbetahat[,j] <- soln$betahat
  QLPqrhat[,j] <- soln$qrhat
  QLPqdifhat[,j] <- soln$qdifhat
  #Estimates from Productivity Differentials
  QLPDSPB[,j] <- soln$DSPB
  QLPDSPC[,j] <- soln$DSPC
  if (j==1){
    LPhat <- soln$LPhat
    LPboot <- soln$LPboot
    #Estimates from Productivity Differentials
    LPPB <- soln$LPPB
    LPPC <- soln$LPPC
    LPPB_boot <- soln$LPPBboot
    LPPC_boot <- soln$LPPCboot
  }
}
filename <- paste("PFQR/DATA/COL/QLP_Environments/QLP_Boot_COL_ISIC", id, ".RData", sep="")
save(tauvec, dZLP, QLPbetahat, LPhat, QLPqrhat, QLPqdifhat, QLPbetaboot, 
  LPboot, QLPqrboot, QLPqdifboot, QLPDSPB, QLPDSPC, 
  QLPDSPB_boot, QLPDSPC_boot, LPPB, LPPC, LPPB_boot, LPPC_boot, file=filename)
#HPC Job Submissions for batches: qsub -t 1:length(industries) DS_Boot_US.job
#Here length(industries)=4


