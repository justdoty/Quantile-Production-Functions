setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
require(prodest)
source('prodestLP.R')

#Simple OLS to obtain estimates

#Load Chilean dataset
chile_panel <- read.csv('chile_panel.csv')
#Choose which industry to select
industries <- c(311, 381, 321, 331)
lm.soln <- list()
for (i in 1:length(industries)){
	chile <- subset(chile_panel, ciiu_3d==industries[i])
	lm.soln[[i]] <-lm(chile$lnva~chile$lnk+chile$lnw+chile$lnb)
}
save(lm.soln, file='OLS_Estimates.RData')
