#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP_Boot.R')
require(stringr)
#Load US dataset
CHL_panel <- read.csv("PFQR/DATA/CHL/CHLdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
CHLdata <- transmute(CHL_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), isic3=isic3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "381", "321", All)
#Vector of quantiles
tau <- seq(5, 95, length.out=19)/100
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tau[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#This gives the "true" estimates using the "true" data
betahat <- array(0, dim=c(dZ, length(industries)))
#This gives the "true" estimates of TFP dispersion using the "true" data
ratiohat <- array(0, dim=c(3, length(industries)))
#Store results for bootstrap replications across quantiles across industries
betaboot <- array(0, dim=c(R, dZ, length(industries)))
#Store results for bootstrap replications across quantiles across industries
ratioboot <- array(0, dim=c(R, 3, length(industries)))

for (isic in 1:length(industries)){
  CHL <- filter(CHLdata, str_detect(isic3, industries[isic]))
  soln <- QLP_Boot(tau=tau, idvar=CHL$id, timevar=CHL$year, Y=CHL$lnva, K=CHL$lnk, L=CHL$lnl, proxy=CHL$lnm, binit=NULL, R=R)
  betahat[,isic] <- soln$betahat
  ratiohat[,isic] <- soln$ratiohat
  betaboot[,,isic] <- soln$betaboot
  ratioboot[,,isic] <- soln$ratiohat
}
filename <- paste("PFQR/DATA/CHL/QLP_Environments/QLP_Boot_CHL_Q", id, ".RData", sep="")
save(betahat, ratiohat, betaboot, ratioboot, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) QLP_Boot_CHL.job
#Here length(tau)=19, recommend splitting into two jobs 1:10 and 11:19




