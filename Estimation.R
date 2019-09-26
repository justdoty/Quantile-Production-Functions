require(prodest)
require(foreign)
setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
chile_panel <- read.csv("chile_panel.csv")
chile <- na.omit(subset(chile_panel, ciiu_3d==381))

#DO NOT USE 311, 331, 322 for ACF

#Use 381, 321 work best for ACF




#Robust check for LP estimation
LP.fit <- prodestLP(Y=chile$lnva, fX=cbind(chile$lnb, chile$lnw), sX=chile$lnk, pX=chile$proxy_m, idvar=chile$id, timevar=chile$year, opt='solnp')
print(LP.fit)
ACF.fit <- prodestACF(Y=chile$lnva, fX=cbind(chile$lnb, chile$lnw), sX=chile$lnk, pX=chile$proxy_e, idvar=chile$id, timevar=chile$year, opt='solnp')
print(ACF.fit)



