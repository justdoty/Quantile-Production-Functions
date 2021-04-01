#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QACF.R')
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
QACFbetahat <- array(0, dim=c(dZ, length(split)))
ACFhat <- array(0, dim=c(dZ, length(split)))
for (t in 1:length(split)){
  COL <- filter(COLdata, year %in% split[[t]])
  soln <- QACF(tau=tau_t, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ, binit=NULL)
  QACFbetahat[,t] <- soln$betahat
  ACFhat[,t] <- soln$ACFhat
}
filename <- paste("PFQR/DATA/COL/QACF/QACF_Environments/QACFT_COL_Q", id, ".RData", sep="")
save(QACFbetahat, ACFhat, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau_t) myjob.job
#Here length(tau_t)=4




