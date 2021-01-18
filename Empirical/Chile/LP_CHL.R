require(stringr)
source('PFQR/FUN/LP.R')

#This function estimates the LP production function batching by industry
#Function calls prodestLP to perform estimation

#Load US dataset
CHL_panel <- read.csv("PFQR/DATA/CHL/CHLdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
CHLdata <- transmute(CHL_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), isic3=isic3)
#Choose which industry to select
All <- "^3"
industries <- c("311", "381", "321", All)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]
#The number of bootstrap replications to be used in LP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Vector of quantiles of TFP
tfptau <- c(0.1, 0.25, 0.5, 0.9)

CHL <- filter(CHLdata, str_detect(isic3, industries))
soln <- LP(idvar=CHL$id, timevar=CHL$year, Y=CHL$lnva, K=CHL$lnk, L=CHL$lnl, proxy=CHL$lnm, dZ=dZ,  binit=NULL, R=R, tfptau=tfptau)
betahat <- soln$betahat
QTFPhat <- soln$QTFPhat
betaboot <- soln$betaboot
QTFPboot <- soln$QTFPboot

filename <- paste("PFQR/DATA/CHL/LP_Environments/LP_CHL_ISIC_", id, ".RData", sep="")
save(betahat, QTFPhat, betaboot, QTFPboot, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4






