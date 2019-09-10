require(foreign)
require(dplyr)
require(prodest)
setwd('/Users/justindoty/Documents/Reading/Econometrics/Microeconometrics/Production/Chile_Data/')
chile <- read.dta('mergecap.dta')

chile <- subset(chile, ciiu_3d==381, select=c(id, year, enter, exit, unsklab, sklab, avgemps, rva, tnk80, ciiu_3d, ciiu3d_descr, 
	relecv, ciiu, tnkall, rmatls))
#relec is ambigious, did LP mean relecv?
#rmatls might be material proxy
names(chile) <- c('ppn', 'year', 'enter','exit', 'blue', 'white', 'L', 'va', 'tnk80', 'ciiu_3d', 'ciiu3d_descr', 'ew1', 'isic', 'tnkall', 'm')

#Drop observations where some of the variables have unreasonable values
chile <- subset(chile, (tnkall>1&!is.na(tnkall)) & va>1 & L>1 & blue<L & blue>0 & white>0 & ew1>0 & m>0)
#Taking logarithms of all the variables in order to do the production function estimation

chile$white <- log(chile$white)
chile$blue <- log(chile$blue)
chile$L <- log(chile$L)
chile$tnkall <- log(chile$tnkall)
chile$va <- log(chile$va)
chile$ew1 <- log(chile$ew1)
chile$m <- log(chile$m)

names(chile) <- c('ppn', 'year', 'enter','exit', 'lnb', 'lnw', 'lnl', 'lnva', 'tnk80', 'ciiu_3d', 'ciiu3d_descr', 'proxy_e', 'isic', 'lnk','proxy_m')

#Deal with plants that enter and exit and drop years

chile <- chile[order(chile$ppn, chile$year), ]

chile <- group_by(chile, ppn) %>% mutate(appear=ifelse(row_number()==1&enter[1]!=1|(row_number()!=1&year[row_number()]-lag(year)!=1&enter[row_number()]!=1),1,0), disappear=ifelse(row_number()==length(ppn)&exit[length(exit)]!=1|(row_number()!=length(ppn)&lead(year)-year[row_number()]!=1&exit[row_number()]!=1),1,0))

chile$appear[chile$year==1979] <- "First Enter"
chile$disappear[chile$year==1996] <- "Last Exit"

chile <- subset(chile, year>1986)

chile <- group_by(chile, ppn) %>% mutate(amil1=max(appear), amil2=max(disappear))
chile <- subset(data.frame(chile), amil1!=1&amil2!=1)



#Robust check for LP estimation
LP.fit <- prodestLP(chile$lnva, fX=cbind(chile$lnb, chile$lnw), sX=chile$lnk, pX=chile$proxy_m, idvar=chile$ppn, timevar=chile$year)
print(LP.fit)
ACF.fit <- prodestACF(chile$lnva, fX=chile$L, sX=chile$lnk, pX=chile$proxy_e, idvar=chile$ppn, timevar=chile$year)
print(ACF.fit)


