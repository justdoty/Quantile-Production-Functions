#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('QLP.R')
require(stringr)
#Load COL dataset
COLdata <- read.csv("COLdata.csv")
#Vector of quantiles
tau_t <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T <- 1
time <- unique(sort(COLdata$year))
split <- split(time, ceiling(seq_along(time)/T))
id <- as.numeric(commandArgs(TRUE)[1])
tau_t <- tau_t[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 10
#The number of parameters being estimated
dZ <- 2
#Bandwidth choice: user specified (for now)
h <- 1e-6
#Store results for bootstrap replications across quantiles across industries
results_T <- array(0, dim=c(R, dZ, length(split)))
#This gives the "true" estimates using the "true" data
true.beta_T <- array(0, dim=c(dZ, length(split)))
for (t in 1:length(split)){
  COL <- filter(COLdata, year %in% split[[t]])
  soln <- tryCatch(QLP(tau=tau_t, va=COL$lnva, state=COL$lnk, free=COL$lnl, proxy=COL$lnm, id=COL$id, time=COL$year, h=h, b.init=NULL, R=R))
  results_T[,,t] <- soln[[1]]
  true.beta_T[,t] <- soln[[2]]
}
filename <- paste("QLPT_COL_Q", id, ".RData", sep="")
save(results_T, true.beta_T, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau_t) myjob.job
#Here length(tau_t)=5




