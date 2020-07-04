# setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
require(prodest)
source('prodestLP.R')

#This function estimates the LP production function batching by industry
#Function calls prodestLP to perform estimation

#Load Chilean dataset
chile_panel <- read.csv('chile_panel.csv')
#Choose which industry to select
industries <- c(311, 381, 321, 331)
id <- as.numeric(commandArgs(TRUE)[1])
industries <- industries[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 3
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ))

chile <- subset(chile_panel, ciiu_3d==industries)
soln <- tryCatch(LP(Y=chile$lnva, sX=chile$lnk, fX=cbind(chile$lnw, chile$lnb), pX=chile$proxy_e, idvar=chile$id, timevar=chile$year, theta0=NULL, R=R))
results <- soln[[1]]
true.beta.LP <- soln[[2]]

filename <- paste("LP_Estimation", id, ".RData", sep="")
save(results, true.beta.LP, file=filename)







