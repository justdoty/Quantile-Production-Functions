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

for (naics in 1:length(industries)){
  US <- filter(USdata, str_detect(naics2, industries[naics]))
  soln <- QLP_Boot(tau=tau, idvar=US$id, timevar=US$year, Y=US$lnva, K=US$lnk, L=US$lnl, proxy=US$lnm, dZ=dZ, binit=NULL, R=R)
  #Bootstrapped Estimates
  QLPbetaboot[,,naics] <- soln$betaboot
  LPboot[,,naics] <- soln$LPboot
  QLPqrboot[,,naics] <- soln$qrboot
  QLPqdifboot[,,naics] <- soln$qdifboot
  #Estimataes from Unsampled data
  QLPbetahat[,naics] <- soln$betahat
  LPhat[,naics] <- soln$LPhat
  QLPqrhat[,naics] <- soln$qrhat
  QLPqdifhat[,naics] <- soln$qdifhat
  QLPTFPhat[[naics]] <- soln$QTFPhat
  LPTFPhat[[naics]] <- soln$TFPhat
  LPomegahat[[naics]] <- soln$omegahat
  LPexpost[[naics]] <- soln$expost
}
filename <- paste("PFQR/DATA/US/QLP_Environments/QLP_Boot_US_Q", id, ".RData", sep="")
save(tauvec, dZ, QLPbetahat, LPhat, QLPqrhat, QLPqdifhat, QLPTFPhat, LPTFPhat, LPomegahat, QLPbetaboot, LPboot, QLPqrboot, QLPqdifboot, LPexpost=LPexpost, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) QLP__Boot_US.job
#Here length(tau)=19, recommend splitting into two jobs 1:10 and 11:19


