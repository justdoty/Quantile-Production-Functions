#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP_Boot.R')
require(stringr)
#Load US dataset
US_panel <- read.csv("PFQR/DATA/US/USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=as.character(naics3), naics2=str_extract(as.character(naics3), "^.{2}"))
#Choose which industry to select
All <- "^3"
industries <- c("31", "32", "33", All)
#Vector of quantiles of firm-size
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

for (naics in 1:length(industries)){
  US <- filter(USdata, str_detect(naics2, industries[naics]))
  soln <- QLP_Boot(tau=tau, idvar=US$id, timevar=US$year, Y=US$lnva, K=US$lnk, L=US$lnl, proxy=US$lnm, dZ=dZ, binit=NULL, R=R, tfptau=tfptau)
  #Bootstrapped Estimates
  betaboot[,,naics] <- soln$betaboot
  qrboot[,,naics] <- soln$qrboot
  qdifboot[,,naics] <- soln$qdifboot
  QTFPboot[,,naics] <- soln$QTFPboot
  #Estimataes from Unsampled data
  betahat[,naics] <- soln$betahat
  qrhat[,naics] <- soln$qrhat
  qdifhat[,naics] <- soln$qdifhat
  QTFPhat[,naics] <- soln$QTFPhat
}
filename <- paste("PFQR/DATA/US/QLP_Environments/QLP_Boot_US_Q", id, ".RData", sep="")
save(tauvec, tfptau, dZ, betahat, qrhat, qdifhat, QTFPhat, betaboot, qrboot, qdifboot, QTFPboot, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) QLP__Boot_US.job
#Here length(tau)=19, recommend splitting into two jobs 1:10 and 11:19


