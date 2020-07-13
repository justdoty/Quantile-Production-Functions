require(stringr)
require(quantreg)
require(dplyr)

#Simple OLS and QR to obtain estimates

#Load US dataset
COLdata <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/COLdata.csv")
#Choose which industry to select
All <- "^3"
industries <- c("311", "322", "381", All)
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
dZ <- 2
lm.soln <- array(0, dim=c(length(industries), dZ))
qr.soln <- array(0, dim=c(length(tau), dZ, length(industries)))
for (isic in 1:length(industries)){
	COL <- filter(COLdata, str_detect(isic3, industries[isic]))
	lm.soln[isic,] <-lm(COL$lnva~COL$lnk+COL$lnl)$coef[-1]
	for (q in 1:length(tau)){
		qr.soln[,,isic][q,] <- rq(COL$lnva~COL$lnk+COL$lnl, tau=tau[q])$coef[2:3]
	}
}
save(lm.soln, qr.soln, file='/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/COL_QR_OLS.Rdata')
