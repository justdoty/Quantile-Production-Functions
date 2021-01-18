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
#The number of bootstrap replications to be used in LP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Vector of quantiles of TFP
tfptau <- c(0.1, 0.25, 0.5, 0.9)


COL <- filter(COLdata, str_detect(isic3, industries))
soln <- LP(idvar=COL$id, timevar=COL$year, Y=COL$lnva, K=COL$lnk, L=COL$lnl, proxy=COL$lnm, dZ=dZ,  binit=NULL, R=R, tfptau=tfptau)
betahat <- soln$betahat
QTFPhat <- soln$QTFPhat
betaboot <- soln$betaboot
QTFPboot <- soln$QTFPboot

filename <- paste("PFQR/DATA/COL/LP_Environments/LP_COL_ISIC_", id, ".RData", sep="")
save(betahat, QTFPhat, betaboot, QTFPboot, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4






