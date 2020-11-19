require(stringr)
source('PFQR/FUN/LP.R')

#This function estimates the LP production function batching by industry
#Function calls prodestLP to perform estimation

#Load US dataset
COLdata <- read.csv("PFQR/DATA/COL/COLdata.csv")
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ))

COL <- filter(COLdata, str_detect(isic3, industries))
soln <- LP(idvar=CHL$id, timevar=CHL$year, Y=CHL$lnva, K=CHL$lnk, L=CHL$lnl, proxy=CHL$lnm,  binit=NULL, R=R)
results <- soln[[1]]
true.beta.LP <- soln[[2]]

filename <- paste("PFQR/DATA/COL/LP_COL_ISIC_", id, ".RData", sep="")
save(results, true.beta.LP, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4






