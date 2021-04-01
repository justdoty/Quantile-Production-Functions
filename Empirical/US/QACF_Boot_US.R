#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QACF_Boot.R')
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

for (naics in 1:length(industries)){
  US <- filter(USdata, str_detect(naics2, industries[naics]))
  soln <- QACF_Boot(tau=tau, idvar=US$id, timevar=US$year, Y=US$lnva, K=US$lnk, L=US$lnl, proxy=US$lnm, dZ=dZ, binit=NULL, R=R)
  #Bootstrapped Estimates
  QACFbetaboot[,,naics] <- soln$betaboot
  ACFboot[,,naics] <- soln$ACFboot
  QACFqrboot[,,naics] <- soln$qrboot
  QACFqdifboot[,,naics] <- soln$qdifboot
  #Estimataes from Unsampled data
  QACFbetahat[,naics] <- soln$betahat
  ACFhat[,naics] <- soln$ACFhat
  QACFqrhat[,naics] <- soln$qrhat
  QACFqdifhat[,naics] <- soln$qdifhat
  QACFTFPhat[[naics]] <- soln$QTFPhat
  ACFTFPhat[[naics]] <- soln$TFPhat
  ACFomegahat[[naics]] <- soln$omegahat
  ACFexpost[[naics]] <- soln$expost
}
filename <- paste("PFQR/DATA/US/QACF/QACF_Environments/QACF_Boot_US_Q", id, ".RData", sep="")
save(tauvec, dZ, QACFbetahat, ACFhat, QACFqrhat, QACFqdifhat, QACFTFPhat, ACFTFPhat, ACFomegahat, QACFbetaboot, ACFboot, QACFqrboot, QACFqdifboot, ACFexpost, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) QACF__Boot_US.job
#Here length(tau)=19, recommend splitting into two jobs 1:10 and 11:19


