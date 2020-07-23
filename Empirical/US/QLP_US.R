#This the file that calls the functions to estimate the production function for the US data
#Runs batches over quantiles for each industry, note that this cannot be performed on personal CPU without
#additional adjustments
source('QLP.R')
require(stringr)
#Load US dataset
US_panel <- read.csv("USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=as.character(naics3), naics2=str_extract(as.character(naics3), "^.{2}"))
#Choose which industry to select
All <- "^3"
industries <- c("31", "32", "33", All)
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

for (naics in 1:length(industries)){
  US <- filter(USdata, str_detect(naics2, industries[naics]))
  soln <- tryCatch(QLP(tau=tau, va=US$lnva, state=US$lnk, free=US$lnl, proxy=US$lnm, id=US$id, time=US$year, h=h, b.init=NULL, R=R))
  results[,,naics] <- soln[[1]]
  true.beta[,naics] <- soln[[2]]
}
filename <- paste("QLP_US_Q", id, ".RData", sep="")
save(results, true.beta, file=filename)


#HPC Job Submissions for batches: qsub -t 1:length(tau) myjob.job
#Here length(tau)=10




