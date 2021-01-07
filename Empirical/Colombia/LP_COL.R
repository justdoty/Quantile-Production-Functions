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


COL <- filter(COLdata, str_detect(isic3, industries))
soln <- LP(idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm,  binit=NULL, R=R)
betahat <- soln$betahat
ratiohat <- soln$ratiohat
betaboot <- soln$betaboot
ratioboot <- soln$ratioboot

filename <- paste("PFQR/DATA/COL/LP_Environments/LP_COL_ISIC_", id, ".RData", sep="")
save(betahat, ratiohat, betaboot, ratioboot, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4






