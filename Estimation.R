# setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
require(prodest)
source('ivqr_gmm.R')

#Load Chilean dataset
chile_panel <- read.csv('chile_panel.csv')
#Choose which industry to select
industries <- c(311, 381, 321, 331)
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
tau_n <- length(tau)
id <- as.numeric(commandArgs(TRUE)[1])
tau <- tau[id]
#The number of bootstrap replications to be used in QLP defined below
R <- 500
#The number of parameters being estimated
dZ <- 3
#Bandwidth choice: user specified (for now)
h <- 1e-6
#Store results for bootstrap replications across quantiles across industries
results <- array(0, dim=c(R, dZ, length(industries)))
true.beta <- array(0, dim=length(industries))


for (ISIC in 1:length(industries)){
  chile <- subset(chile_panel, ciiu_3d==industries[ISIC])
  soln <- tryCatch(QLP(tau=tau, va=chile$lnva, state=chile$lnk, free=cbind(chile$lnw, chile$lnb), proxy=chile$proxy_e, id=chile$id, time=chile$year, h=h, b.init=NULL, R=R))
  results[,,ISIC] <- soln[[1]]
  true.beta[ISIC] <- soln[[2]]
}
true.
filename <- paste("QLP_Estimation", id, ".RData", sep="")
save(results, true.beta, file=filename)







