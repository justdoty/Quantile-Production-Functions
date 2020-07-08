setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US')
require(stringr)
require(quantreg)

#Simple OLS and QR to obtain estimates

#Load US dataset
US_panel <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/USdata.csv")
#Convert 3 digit NAICS code to 2 digit NAICS and take natural logs
USdata <- transmute(US_panel, id=id, year=year, lny=log(Y), lnva=log(VA), lnk=log(K), lnl=log(L), lnm=log(M), naics3=naics3, naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Choose which industry to select
industries <- c(31, 32, 33)
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
dZ <- 2
lm.soln <- array(0, dim=c(length(industries), dZ))
qr.soln <- array(0, dim=c(length(tau), dZ, length(industries)))
for (naics in 1:length(industries)){
	US <- subset(USdata, naics2==industries[naics])
	lm.soln[naics,] <-lm(US$lnva~US$lnk+US$lnl)$coef[-1]
	for (q in 1:length(tau)){
		qr.soln[,,naics][q,] <- rq(US$lnva~US$lnk+US$lnl, tau=tau[q])$coef[2:3]
	}
}
save(lm.soln, qr.soln, file='US_QR_OLS.Rdata')
