require(stringr)
source('PFQR/FUN/LP.R')

#This function estimates the LP production function batching by industry

#Load US dataset
US_panel <- read.csv("PFQR/DATA/US/USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Choose which industry to select
All <- "^3"
industries <- c("31", "32", "33", All)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]
#The number of bootstrap replications to be used in LP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2

US <- filter(USdata, str_detect(naics2, industries))
soln <- LP(idvar=US$id, timevar=US$year, Y=US$lnva, K=US$lnk, L=US$lnl, proxy=US$lnm,  binit=NULL, R=R)
betahat <- soln$betahat
ratiohat <- soln$ratiohat
betaboot <- soln$betaboot
ratioboot <- soln$ratioboot

filename <- paste("PFQR/DATA/US/LP_Environments/LP_US_NAICS_", id, ".RData", sep="")
save(betahat, ratiohat, betaboot, ratioboot, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=4






