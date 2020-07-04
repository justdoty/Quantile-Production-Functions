require(prodest)
require(stringr)
source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/prodestLP.R')

#Simple OLS to obtain estimates

#Load US dataset
US_panel <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Choose which industry to select
industries <- c(31, 32, 33)
lm.soln <- list()
for (i in 1:length(industries)){
	US <- subset(USdata, naics2==industries[naics])
	lm.soln[[i]] <-lm(US$lnva~US$lnk+US$lnl)
}
setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US')
save(lm.soln, file='OLS_Estimates_US.RData')
