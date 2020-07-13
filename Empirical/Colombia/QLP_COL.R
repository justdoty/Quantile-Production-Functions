#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('QLP.R')
require(stringr)
#Load US dataset
COLdata <- read.csv("COLdata.csv")
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tau[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Bandwidth choice: user specified (for now)
h <- 1e-6
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ, length(industries)))
#This gives the "true" estimates using the "true" data
true.beta <- array(0, dim=c(dZ, length(industries)))


for (isic in 1:length(industries)){
  COL <- filter(COLdata, str_detect(isic3, industries[isic]))
  soln <- tryCatch(QLP(tau=tau, va=COL$lnva, state=COL$lnk, free=COL$lnl, proxy=COL$lnm, id=COL$id, time=COL$year, h=h, b.init=NULL, R=R))
  results[,,isic] <- soln[[1]]
  true.beta[,isic] <- soln[[2]]
}
filename <- paste("QLP_COL_Q", tau[id], ".RData", sep="")
save(results, true.beta, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) myjob.job
#Here length(tau)=10




