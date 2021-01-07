#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP.R')
require(stringr)
#Load COL dataset
COLdata <- read.csv("PFQR/DATA/COL/COLdata.csv")
#Vector of quantiles
tau_t <- c(0.1, 0.25, 0.5, 0.9)
T <- 2
time <- unique(sort(COLdata$year))
split <- split(time, ceiling(seq_along(time)/T))
id <- as.numeric(commandArgs(TRUE)[1])
tau_t <- tau_t[id]
#The number of parameters being estimated
dZ <- 2
betahat <- array(0, dim=c(dZ, length(split)))
ratiohat <- array(0, dim=c(3, length(split)))
for (t in 1:length(split)){
  COL <- filter(COLdata, year %in% split[[t]])
  soln <- QLP(tau=tau_t, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, binit=NULL)
  betahat[,t] <- soln$betahat
  ratiohat[,t] <- soln$ratiohat
}
filename <- paste("PFQR/DATA/COL/QLP_Environments/Time_Estimates/QLPT_COL_Q", id, ".RData", sep="")
save(betahat, ratiohat, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau_t) myjob.job
#Here length(tau_t)=5




