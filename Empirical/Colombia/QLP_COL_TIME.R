#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP.R')
require(stringr)
#Load COL dataset
COLdata <- read.csv("PFQR/DATA/COL/COLdata.csv")
#Vector of quantiles
tau_t <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T <- 2
time <- unique(sort(COLdata$year))
split <- split(time, ceiling(seq_along(time)/T))
id <- as.numeric(commandArgs(TRUE)[1])
tau_t <- tau_t[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Store results for bootstrap replications across quantiles across industries
results_T <- array(0, dim=c(R, dZ, length(split)))
#This gives the "true" estimates using the "true" data
true.beta_T <- array(0, dim=c(dZ, length(split)))
for (t in 1:length(split)){
  COL <- filter(COLdata, year %in% split[[t]])
  soln <- QLP(tau=tau_t, idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, binit=NULL, R=R)
  results_T[,,t] <- soln[[1]]
  true.beta_T[,t] <- soln[[2]]
}
filename <- paste("PFQR/DATA/COL/QLPT_COL_Q", id, ".RData", sep="")
save(results_T, true.beta_T, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau_t) myjob.job
#Here length(tau_t)=5




