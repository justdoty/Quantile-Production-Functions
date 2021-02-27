#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP_Boot.R')
require(stringr)
#Load US dataset
CHL_panel <- read.csv("PFQR/DATA/CHL/CHLdata.csv")
#Convert 3 digit isic code to 2 digit isic and take natural logs
CHLdata <- transmute(CHL_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), isic3=isic3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "381", "321", All)
#Vector of quantiles
tauvec <- seq(5, 95, length.out=19)/100
#Vector of quantiles of TFP
tfptau <- seq(5, 95, length.out=19)/100
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tauvec[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Bandwidth specification
h <- NULL
#This gives the estimates using the unsampled data########################
#Elasticities from QLP
betahat <- array(0, dim=c(dZ, length(industries)))
LPhat <- array(0, dim=c(dZ, length(industries)))
#Elasticities from QR
qrhat <- array(0, dim=c(dZ, length(industries)))
#Difference between QLP and QR
qdifhat <- array(0, dim=c(dZ, length(industries)))
#Quantiles of TFP
QTFPhat <- array(0, dim=c(length(tfptau), length(industries)))
#Store results for bootstrap replications##############################
#Bootstrapped Elasticiteis from QLP
betaboot <- array(0, dim=c(R, dZ, length(industries)))
LPboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped Elasticities from QR
qrboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped differences between QLP and QR
qdifboot <- array(0, dim=c(R, dZ, length(industries)))
#Bootstrapped quantiles of TFP
QTFPboot <- array(0, dim=c(R, length(tfptau), length(industries)))

for (isic in 1:length(industries)){
  CHL <- filter(CHLdata, str_detect(isic3, industries[isic]))
  soln <- QLP_Boot(tau=tau, h=h, idvar=CHL$id, timevar=CHL$year, Y=CHL$lnva, K=CHL$lnk, L=CHL$lnl, proxy=CHL$lnm, dZ=dZ, binit=NULL, R=R, tfptau=tfptau)
  #Bootstrapped Estimates
  betaboot[,,isic] <- soln$betaboot
  LPboot[,,isic] <- soln$LPboot
  qrboot[,,isic] <- soln$qrboot
  qdifboot[,,isic] <- soln$qdifboot
  QTFPboot[,,isic] <- soln$QTFPboot
  #Estimataes from Unsampled data
  betahat[,isic] <- soln$betahat
  LPhat[,isic] <- soln$LPhat
  qrhat[,isic] <- soln$qrhat
  qdifhat[,isic] <- soln$qdifhat
  QTFPhat[,isic] <- soln$QTFPhat
}
filename <- paste("PFQR/DATA/CHL/QLP_Environments/QLP_Boot_CHL_Q", id, ".RData", sep="")
save(tauvec, tfptau, dZ, betahat, LPhat, qrhat, qdifhat, QTFPhat, betaboot, LPboot, qrboot, qdifboot, QTFPboot, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) QLP_Boot_CHL.job
#Here length(tau)=19, recommend splitting into two jobs 1:10 and 11:19




