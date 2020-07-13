require(prodest)
require(stringr)
source('prodestLP.R')

#This function estimates the LP production function batching by industry
#Function calls prodestLP to perform estimation

#Load US dataset
COLdata <- read.csv("COLdata.csv")
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
soln <- tryCatch(LP(Y=COL$lnva, sX=COL$lnk, fX=COL$lnl, pX=COL$lnm, idvar=COL$id, timevar=COL$year, theta0=NULL, R=R))
results <- soln[[1]]
true.beta.LP <- soln[[2]]

filename <- paste("LP_COL_ISIC_", id, ".RData", sep="")
save(results, true.beta.LP, file=filename)

#HPC Job Submissions for batches: qsub -t 1:length(industries) myjob.job
#Here length(industries)=3






