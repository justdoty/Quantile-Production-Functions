#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('QLP.R')
require(stringr)
#Load US dataset
US_panel <- read.csv("USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=as.character(naics3), naics2=str_extract(as.character(naics3), "^.{2}"))
#Vector of quantiles
tau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
allyears <- unique(USdata$year)
panel_T <- length(allyears)
split_T <- 5
yearsplit <- split(allyears, ceiling(seq_along(allyears)/split_T))
id <- as.numeric(commandArgs(TRUE)[1])
ysplit <- yearsplit[[id]]
#The number of bootstrap replications to be used in QLP defined below
R <- 1
#The number of parameters being estimated
dZ <- 2
#Bandwidth choice: user specified (for now)
h <- 1e-6
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ, length(tau)))
#This gives the "true" estimates using the "true" data
true.beta <- array(0, dim=c(dZ, length(tau)))


for (q in 1:length(tau)){
  US <- filter(USdata, year %in% ysplit)
  soln <- tryCatch(QLP(tau=tau[q], va=US$lnva, state=US$lnk, free=US$lnl, proxy=US$lnm, id=US$id, time=US$year, h=h, b.init=NULL, R=R))
  results_T[,,q] <- soln[[1]]
  true.beta_T[,q] <- soln[[2]]
}
filename <- paste("QLP_US_T", id, ".RData", sep="")
save(results_T, true.beta_T, file=filename)


#HPC Job Submissions for batches: qsub -t 1:T myjob.job
#Here T=length(allyears)/split_T=10




