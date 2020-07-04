require(prodest)
require(stringr)
source('prodestLP.R')

#This function estimates the LP production function batching by industry
#Function calls prodestLP to perform estimation

#Load US dataset
US_panel <- read.csv("USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Choose which industry to select
industries <- c(31, 32, 33)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 2
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ))

US <- subset(USdata, naics2==industries)
soln <- tryCatch(LP(Y=US$lnva, sX=US$lnk, fX=US$lnl, pX=US$lnm, idvar=US$id, timevar=US$year, theta0=NULL, R=R))
#This gives a list of the bootstrapped estimates across quantiles across industries
results <- soln[[1]]
#This gives the "true" estimates using the "true" data
true.beta.LP <- soln[[2]]

filename <- paste("US_LP_Estimation", id, ".RData", sep="")
save(results, true.beta.LP, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=3






