#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('PFQR/FUN/QLP.R')
require(stringr)
#Load US dataset
US_panel <- read.csv("PFQR/DATA/US/USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=as.character(naics3), naics2=str_extract(as.character(naics3), "^.{2}"))
#Vector of quantiles
tau_t <- c(0.1, 0.25, 0.5, 0.9)
T <- 5
time <- unique(sort(USdata$year))
split <- split(time, ceiling(seq_along(time)/T))
id <- as.numeric(commandArgs(TRUE)[1])
tau_t <- tau_t[id]
#The number of parameters being estimated
dZ <- 2
h <- NULL
betahat <- array(0, dim=c(dZ, length(split)))
LPhat <- array(0, dim=c(dZ, length(split)))
for (t in 1:length(split)){
  US <- filter(USdata, year %in% split[[t]])
  soln <- QLP(tau=tau_t, h=h, idvar=US$id, timevar=US$year, Y=US$lnva, K=US$lnk, L=US$lnl, proxy=US$lnm, dZ=dZ, binit=NULL)
  betahat[,t] <- soln$betahat
  LPhat[,t] <- soln$LPhat
}
filename <- paste("PFQR/DATA/US/QLP_Environments/Time_Estimates/QLPT_US_Q", id, ".RData", sep="")
save(betahat, LPhat, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau_t) QLP_US_TIME.job
#Here length(tau_t)=4




